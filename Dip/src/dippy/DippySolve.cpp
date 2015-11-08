#include "Python.h"

#include "UtilParameters.h"
//===========================================================================//
#include "DippyDecompApp.h"
#include "DippyDecompAlgo.h"
#include "DippyPythonUtils.h"
//===========================================================================//
#include "AlpsDecompModel.h"
//===========================================================================//
#include "DecompAlgoC.h"
#include "DecompAlgoD.h"
#include "DecompAlgoPC.h"
#include "DecompAlgoRC.h"
#include "UtilMacros.h"
//===========================================================================//

#include <assert.h>
#ifdef _DEBUG
#ifdef USEVLD
#include <vld.h>
#endif
#endif

#ifdef _WIN32
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT extern "C"
#endif

#ifndef DECOMP_INF_DEFINED
#define DECOMP_INF_DEFINED
double DecompInf = COIN_DBL_MAX;
#else
extern double DecompInf;
#endif

DLLEXPORT PyObject* Solve(PyObject* self, PyObject* args)
{
   PyObject* pProb;
   PyObject* pParamDict;

   if (!PyArg_ParseTuple(args, "OO", &pProb, &pParamDict)) {
      return NULL;
   }

   try {
      // create the utility class for parsing parameters
      UtilParameters utilParam;
      // By default Dippy enforces branching in the master problem
      utilParam.Add("DECOMP", "BranchEnforceInMaster", "1");
      utilParam.Add("DECOMP", "BranchEnforceInSubProb", "0");
      // loop through paramDict, add to utilParam
      PyObject* pKey, *pValue;
      Py_ssize_t pos = 0;

      while (PyDict_Next(pParamDict, &pos, &pKey, &pValue)) {
         // key is a 2-tuple (section, name), both strings
         // value is a string
         // TODO: better error reporting
         const char* section = NULL;
         PyObject* pSection = PyTuple_GetItem(pKey, 0);

         if (pSection != Py_None) {
            section = PyString_AsString(PyTuple_GetItem(pKey, 0));
         }

         const char* name = PyString_AsString(PyTuple_GetItem(pKey, 1));
         const char* value = PyString_AsString(pValue);
         utilParam.Add(section, name, value);
      }

      bool doCut = utilParam.GetSetting("doCut", false);
      bool doPriceCut = utilParam.GetSetting("doPriceCut", false);
      bool doRelaxCut = utilParam.GetSetting("doRelaxCut", false);
      // create the user application (a DecompApp)
      DippyDecompApp sip(utilParam, pProb);
      // create the decomp algo
      DecompAlgo* algo = NULL;

      if (doPriceCut) {
         algo = new DippyAlgoPC(&sip, utilParam, pProb);
      } else if (doCut) {
         algo = new DippyAlgoC(&sip, utilParam, pProb);
      } else if (doRelaxCut) {
         algo = new DippyAlgoRC(&sip, utilParam, pProb);
      }

      // default
      if (algo == NULL) {
         algo = new DippyAlgoC(&sip, utilParam, pProb);
      }

      AlpsDecompModel alpsModel(utilParam, algo);
      alpsModel.solve();
      // TODO: Python exception needs to be set here or higher
      int status = alpsModel.getSolStatus();
      PyObject* pStatus;
      PyObject* pMessage = Py_None;
      /**
        LpStatusOptimal     “Optimal”      1
        LpStatusNotSolved   “Not Solved”   0
        LpStatusInfeasible  “Infeasible”  -1
        LpStatusUnbounded   “Unbounded”   -2
        LpStatusUndefined   “Undefined”   -3
      */

      switch (status) {
      case AlpsExitStatusOptimal:
         pStatus = PyInt_FromLong(1);
         break;

      case AlpsExitStatusTimeLimit:
         pStatus = PyInt_FromLong(0);
         pMessage = PyString_FromString("Reached time limit");
         break;

      case AlpsExitStatusNodeLimit:
         pStatus = PyInt_FromLong(0);
         pMessage = PyString_FromString("Reached node limit");
         break;

      case AlpsExitStatusSolLimit:
         pStatus = PyInt_FromLong(0);
         pMessage = PyString_FromString("Reached sol limit");
         break;

      case AlpsExitStatusInfeasible:
         pStatus = PyInt_FromLong(-1);
         break;

      case AlpsExitStatusNoMemory:
         throw UtilException("Out of memory", "Solve", "DippySolve");

      case AlpsExitStatusFailed:
         throw UtilException("Solve failed", "Solve", "DippySolve");

      case AlpsExitStatusUnbounded:
         pStatus = PyInt_FromLong(-2);
         break;

      case AlpsExitStatusFeasible:
         throw UtilException("Feasible but not optimal", "Solve", "DippySolve");

      default:
         throw UtilException("Unknown solution status", "Solve", "DippySolve");
      }

      const DecompSolution* solution = alpsModel.getBestSolution();
      //        cout << "Optimal Solution" << endl;
      //        solution->print();
      PyObject* pSolution = Py_None;

      if (solution != NULL) {
         const double* values = solution->getValues();
         pSolution = pyTupleList_FromDoubleArray(values, sip.m_colList);
      }

      PyObject* pDuals = Py_None;

      if (doCut) {
         DecompAlgoC* algoC = dynamic_cast<DecompAlgoC*>(algo);
         OsiSolverInterface* masterOSI = algoC->getMasterOSI();
         const double* duals = masterOSI->getRowPrice();

         if (duals != NULL) {
            pDuals = pyTupleList_FromDoubleArray(duals, sip.m_rowList);
         }
      }

      delete algo;
      // return solution as list
      PyObject* pOutput = PyTuple_New(4);
      PyTuple_SetItem(pOutput, 0, pStatus);
      PyTuple_SetItem(pOutput, 1, pMessage);
      PyTuple_SetItem(pOutput, 2, pSolution);
      PyTuple_SetItem(pOutput, 3, pDuals);
      Py_INCREF(pOutput);
      return pOutput;
   } catch (CoinError& ex) {
      cerr << "COIN Exception [ " << ex.message() << " ]"
           << " at " << ex.fileName()  << ":L" << ex.lineNumber()
           << " in " << ex.className() << "::" << ex.methodName() << endl;
      return NULL;
   }

   Py_INCREF(Py_None);
   return Py_None;
}
