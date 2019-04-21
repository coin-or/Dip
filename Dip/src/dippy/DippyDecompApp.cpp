#include "Python.h"

#include "DecompVar.h"
#include "DecompAlgoC.h"
#include "DippyDecompApp.h"
#include "DippyDecompCut.h"
#include "DippyPythonUtils.h"

#include "UtilMacros.h"
#include "UtilMacrosDecomp.h"

//===========================================================================//

/**
 * Called to create the core and relaxation models
 */
void DippyDecompApp::createModels()
{
   int i, len;
   string name;
   // create the core model
   DecompConstraintSet* modelCore = new DecompConstraintSet();
   // gets the master problem model
   char masterTuple[] = "getMasterAsTuple";
   PyObject* pMasterAsTuple = PyObject_CallMethod(m_pProb, masterTuple, NULL);

   if (pMasterAsTuple == NULL) {
      throw UtilException("Error calling method prob.getMasterAsTuple()", 
			  "createModels", "DippyDecompApp");
   }

   PyObject* pObjective = PyTuple_GetItem(pMasterAsTuple, 0);
   PyObject* pRowList   = PyTuple_GetItem(pMasterAsTuple, 1);
   PyObject* pColList   = PyTuple_GetItem(pMasterAsTuple, 2);
   m_rowList = pRowList;
   Py_XINCREF(m_rowList);
   m_numCols = PyObject_Length(pColList);
   m_colList = pColList;
   Py_XINCREF(m_colList);
   int numRows = PyObject_Length(pRowList);
   PyObject* pRow, *pRowName, *pRowLb, *pRowUb;
   double lb, ub;

   for (int i = 0; i < numRows; i++) {
      pRow = PyList_GetItem(pRowList, i);
      char getName[] = "getName";
      pRowName = PyObject_CallMethod(pRow, getName, NULL);

      if (pRowName == NULL) {
         throw UtilException("Error calling method row.getName()", 
			     "createModels", "DippyDecompApp");
      }
      char getLb[] = "getLb";
      pRowLb   = PyObject_CallMethod(pRow, getLb, NULL);

      if (pRowLb == NULL) {
         throw UtilException("Error calling method row.getLb()", 
			     "createModels", "DippyDecompApp");
      }
      char getUb[] = "getUb";
      pRowUb   = PyObject_CallMethod(pRow, getUb, NULL);

      if (pRowUb == NULL) {
         throw UtilException("Error calling method row.getUb()", 
			     "createModels", "DippyDecompApp");
      }

      name = PyBytes_AsString(PyUnicode_AsEncodedString(pRowName, "UTF-8", "strict"));

      if (pRowLb == Py_None) {
         lb = -m_infinity;
      } else {
         lb = PyFloat_AsDouble(pRowLb);
      }

      if (pRowUb == Py_None) {
         ub = m_infinity;
      } else {
         ub = PyFloat_AsDouble(pRowUb);
      }

      modelCore->rowNames.push_back(name);
      modelCore->rowLB.push_back(lb);
      modelCore->rowUB.push_back(ub);
      m_rowIndices[pRow] = i; // Don't need to increase reference count here as m_rowList
      // references pRow
   }

   PyObject* pCol, *pColName, *pColLb, *pColUb, *pIsInt;

   for (int j = 0; j < m_numCols; j++) {
      pCol = PyList_GetItem(pColList, j);
      char getName[] = "getName";
      pColName = PyObject_CallMethod(pCol, getName, NULL);

      if (pColName == NULL) {
         throw UtilException("Error calling method col.getName()", 
			     "createModels", "DippyDecompApp");
      }
      char getLb[] = "getLb";
      pColLb   = PyObject_CallMethod(pCol, getLb, NULL);

      if (pColLb == NULL) {
         throw UtilException("Error calling method col.getLb()", 
			     "createModels", "DippyDecompApp");
      }
      char getUb[] = "getUb";
      pColUb   = PyObject_CallMethod(pCol, getUb, NULL);

      if (pColUb == NULL) {
         throw UtilException("Error calling method col.getUb()", 
			     "createModels", "DippyDecompApp");
      }
      char isInteger[] = "isInteger";
      pIsInt   = PyObject_CallMethod(pCol, isInteger, NULL);

      if (pIsInt == NULL) {
         throw UtilException("Error calling method col.isInteger()", 
			     "createModels", "DippyDecompApp");
      }

      name = PyBytes_AsString(PyUnicode_AsEncodedString(pColName, "UTF-8", "strict"));

      if (pColLb == Py_None) {
         lb = -m_infinity;
      } else {
         lb = PyFloat_AsDouble(pColLb);
      }

      if (pColUb == Py_None) {
         ub = m_infinity;
      } else {
         ub = PyFloat_AsDouble(pColUb);
      }

      modelCore->colNames.push_back(name);
      modelCore->colLB.push_back(lb);
      modelCore->colUB.push_back(ub);

      if (PyObject_IsTrue(pIsInt)) {
         modelCore->integerVars.push_back(j);
      }

      m_colIndices[pCol] = j; // Don't need to increase reference count here 
      // as m_rowList references pCol
   }

   // set objective coefficients
   double* obj = new double[m_numCols];
   UtilFillN(obj, m_numCols, 0.0);
   PyObject* pObjKeys = PyDict_Keys(pObjective);
   PyObject* pCoeff;

   for (int j = 0; j < PyObject_Length(pObjKeys); j++) {
      pCol = PyList_GetItem(pObjKeys, j);
      pCoeff = PyDict_GetItem(pObjective, pCol);
      obj[m_colIndices[pCol]] = PyFloat_AsDouble(pCoeff);
   }

   setModelObjective(obj, m_numCols);
   // set constraint matrix
   modelCore->M = pyConstraints_AsPackedMatrix(pRowList, m_rowIndices, 
					       m_colIndices);
   modelCore->M->setDimensions(modelCore->rowLB.size(), 
			       modelCore->colLB.size());
   // subproblems
   char getRelaxsAsDict[] = "getRelaxsAsDict";
   PyObject* pRelaxedDict = PyObject_CallMethod(m_pProb, getRelaxsAsDict, NULL);

   if (pRelaxedDict == NULL) {
      throw UtilException("Error calling method prob.getRelaxsAsDict()", 
			  "createModels", "DippyDecompApp");
   }

   int* masterOnly = new int[m_numCols];

   if (!masterOnly) {
      throw UtilExceptionMemory("createModels", "DecompApp");
   }

   UtilFillN(masterOnly, m_numCols, 1);

   int nRelaxed = 0;

   if (pRelaxedDict != Py_None) {
      nRelaxed = PyObject_Length(pRelaxedDict);
   }

   // we have a list of relaxations
   m_relaxedKeys = PyDict_Keys(pRelaxedDict);
   Py_XINCREF(m_relaxedKeys);
   PyObject* pKey, *pRelax;

   for (int p = 0; p < nRelaxed; p++) {
      DecompConstraintSet* modelRelax = new DecompConstraintSet();
      // each relaxation is a LpProblem
      pKey = PyList_GetItem(m_relaxedKeys, p);
      pRelax = PyDict_GetItem(pRelaxedDict, pKey);
      m_relaxIndices[pKey] = p; // Don't need to increase reference count here 
      //as m_relaxedKey references pKey
      char getRelaxAsTuple[] = "getRelaxAsTuple";
      char O[] = "O";
      PyObject* pRelaxAsTuple = PyObject_CallMethod(m_pProb,getRelaxAsTuple,
                                                    O, pRelax);

      if (pRelaxAsTuple == NULL) {
         throw UtilException("Error calling method prob.getRelaxAsTuple()", 
			     "createModels", "DippyDecompApp");
      }

      // row names
      pRowList = PyTuple_GetItem(pRelaxAsTuple, 0);
      pColList = PyTuple_GetItem(pRelaxAsTuple, 1);
      numRows = PyObject_Length(pRowList);
      map<PyObject*, int> relaxRowIndices;

      for (int i = 0; i < numRows; i++) {
         pRow = PyList_GetItem(pRowList, i);
	 char getName[] = "getName";
         pRowName = PyObject_CallMethod(pRow, getName, NULL);

         if (pRowName == NULL) {
            throw UtilException("Error calling method row.getName()", 
				"createModels", "DippyDecompApp");
         }
	 char getLb[] = "getLb";
         pRowLb   = PyObject_CallMethod(pRow, getLb, NULL);

         if (pRowLb == NULL) {
            throw UtilException("Error calling method row.getLb()", 
				"createModels", "DippyDecompApp");
         }
	 char getUb[] = "getUb";
         pRowUb   = PyObject_CallMethod(pRow, getUb, NULL);

         if (pRowUb == NULL) {
            throw UtilException("Error calling method row.getUb()", 
				"createModels", "DippyDecompApp");
         }

         name = PyBytes_AsString(PyUnicode_AsEncodedString(pRowName, "UTF-8", "strict"));

         if (pRowLb == Py_None) {
            lb = -m_infinity;
         } else {
            lb = PyFloat_AsDouble(pRowLb);
         }

         if (pRowUb == Py_None) {
            ub = m_infinity;
         } else {
            ub = PyFloat_AsDouble(pRowUb);
         }

         modelRelax->rowNames.push_back(name);
         modelRelax->rowLB.push_back(lb);
         modelRelax->rowUB.push_back(ub);
         relaxRowIndices[pRow] = i;
      }

      // get the constraint matrix for this relaxation
      modelRelax->M = pyConstraints_AsPackedMatrix(pRowList, relaxRowIndices,
						   m_colIndices);  

      // set all cols at their lower bounds
      for (int j = 0; j < modelCore->colLB.size(); j++) {
         modelRelax->colLB.push_back(modelCore->colLB[j]);
         modelRelax->colUB.push_back(modelCore->colLB[j]);
      }

      // get active cols
      int cols = PyObject_Length(pColList);
      int index;

      for (int j = 0; j < cols; j++) {
         pCol = PyList_GetItem(pColList, j);
         index = m_colIndices[pCol];

         if ( (index < 0) || (index >= m_colIndices.size()) ) {
            throw UtilException("Bad index for " + name, "createModels", 
				"DippyDecompApp");
         }

         modelRelax->colUB[index] = modelCore->colUB[index];
         modelRelax->activeColumns.push_back(index);
         masterOnly[index] = 0;
      }

      modelRelax->M->setDimensions(modelRelax->rowLB.size(), 
				   modelRelax->colLB.size());

      // copy integer vars (from master prob)
      for (int j = 0; j < modelCore->integerVars.size(); j++) {
         modelRelax->integerVars.push_back(modelCore->integerVars[j]);
      }

      setModelRelax(modelRelax, "BLOCK", p);
   }

   for (i = 0; i < m_numCols; i++) {
      if (masterOnly[i]){
	 modelCore->masterOnlyCols.push_back(i);
      }
   }

   printf("Num master-only cols: %d\n", modelCore->masterOnlyCols.size());

   // set the core problem
   setModelCore(modelCore, "CORE");
   UTIL_DELARR(masterOnly);
   
   assert(!PyErr_Occurred());
}

/**
 * solveRelaxed callback
 *
 * This is called by DIP. This function interfaces with Python to
 * call the user defined function if it's present
 *
 * We're expected to populate varList (basically a vector) and
 * return the status of the subproblem solver.
 */
DecompSolverStatus DippyDecompApp::solveRelaxed(const int whichBlock,
						const double* redCostX, 
						const double convexDual, 
						DecompVarList& varList)
{
   if (!m_pySolveRelaxed) {
      return DecompSolStatNoSolution;
   }

   PyObject* pRelaxKey = PyList_GetItem(m_relaxedKeys, whichBlock);
   PyObject* pRedCostList = pyTupleList_FromDoubleArray(redCostX, m_colList);
   PyObject* pConvexDual = PyFloat_FromDouble(convexDual);
   // call solveRelaxed on DipProblem

   char solveRelaxed[] = "solveRelaxed";
   char OOd[] = "OOd";
   PyObject* pStatandVarList = PyObject_CallMethod(m_pProb, solveRelaxed, OOd, 
					             pRelaxKey,
					             pRedCostList,
					             pConvexDual);

   Py_DECREF(pRedCostList);
   Py_DECREF(pConvexDual);

   if ( (pStatandVarList == NULL) || (pStatandVarList == Py_None) ){
      throw UtilException("Error calling method prob.solveRelaxed()", "solveRelaxed",
                          "DippyDecompApp");
   }

   // [status, varList] = relaxed_solver(...)
   PyObject * pStatus = PyTuple_GetItem(pStatandVarList, 0);

   int cStatus = PyLong_AsLong(pStatus);

   DecompSolverStatus status = (DecompSolverStatus)cStatus;

   PyObject * pVarList = PyTuple_GetItem(pStatandVarList, 1);

   int nVars = PyObject_Length(pVarList);

   // In the new design, we need to allow the possibility that the user will solve
   // the problem exactly, but not find any solutions with reduced costs zero
   // The below is is commented out and left in the source for posterity
   // tkr 11/11/15
   //if (nVars == 0) {
   //   throw UtilException("Empty variable list", "solveRelaxed", "DippyDecompApp");
   //}

   // solveRelaxed returns 3-tuples (cost, reduced cost, dictionary of (variable, value) pairs)
   // We can use these to construct a C++ DecompVar objects
   double cost, rc;
   PyObject* pTuple, *pDict, *pKeys, *pCol;
   string name;
   double value;

   for (int j = 0; j < nVars; j++) {
      pTuple = PySequence_GetItem(pVarList, j);
      cost   = PyFloat_AsDouble(PyTuple_GetItem(pTuple, 0));
      rc     = PyFloat_AsDouble(PyTuple_GetItem(pTuple, 1));

      pDict = PyTuple_GetItem(pTuple, 2);
      pKeys = PyDict_Keys(pDict);
      vector<int>    varInds;
      vector<double> varVals;

      for (int n = 0; n < PyObject_Length(pDict); n++) {
         pCol  = PyList_GetItem(pKeys, n);
         value = PyFloat_AsDouble(PyDict_GetItem(pDict, pCol));
         varInds.push_back(m_colIndices[pCol]);
         varVals.push_back(value);
      }

      Py_DECREF(pKeys);
      Py_DECREF(pTuple);

      DecompVar* var =  new DecompVar(varInds, varVals, rc, cost);
      var->setBlockId(whichBlock);
      varList.push_back(var);
   }

   Py_DECREF(pStatandVarList);

   assert(!PyErr_Occurred());

   return status;
}

/**
 * APPisUserFeasible callback
 *
 * Called by DIP, we interface with Python
 */
bool DippyDecompApp::APPisUserFeasible(const double* x, const int n_cols, const double tolZero)
{
   assert(n_cols == m_modelCore.getModel()->getColNames().size());
   PyObject* pSolutionList = pyTupleList_FromDoubleArray(x, m_colList);
   PyObject* pTolZero = PyFloat_FromDouble(tolZero);

   if (!m_pyIsUserFeasible) {
      return true;
   }
   char isUserFeasible[] = "isUserFeasible";
   char Od[] = "Od";
   PyObject* pResult = PyObject_CallMethod(m_pProb, isUserFeasible, Od, pSolutionList, pTolZero);

   if (pResult == NULL) {
      throw UtilException("Error calling method prob.isUserFeasible()", "APPisUserFeasible",
                          "DippyDecompApp");
   }

   // This should not happen as having isUserFeasible present but not returning a boolean is
   // not good
   if (pResult == Py_None) {
      // method exists, but not implemented, return true
      return true;
   }

   assert(!PyErr_Occurred());

   return (bool)PyObject_IsTrue(pResult);
}

/**
 * generateCuts callback
 *
 * Called by DIP, we interface with Python
 */

int DippyDecompApp::generateCuts(const double* x, DecompCutList& cutList)
{
   if (!m_pyGenerateCuts) {
      return 0;
   }

   // PyObject *pSolutionList = pyTupleList_FromDoubleArray(x, m_colList);
   // MO (28/2/2012) - Don't need this anymore as solution is contained within node
   PyObject* pPackagedNode = pyTupleList_FromNode(getDecompAlgo(), STAT_FEASIBLE);
   char generateCuts[] = "generateCuts";
   char O[] = "O";
   PyObject* pCutList = PyObject_CallMethod(m_pProb, generateCuts, O, pPackagedNode);

   if (pCutList == NULL) {
      throw UtilException("Error calling method prob.generateCuts()", "generateCuts",
                          "DippyDecompApp");
   }

   // This should never happen, pyGenerateCuts should be set to false in dippy.py
   if (pCutList == Py_None)
      // method exists, but is not implemented, return 0
   {
      return 0;
   }

   // generateCuts returns constraints, i.e., dictionary of (variable, value) pairs also with name, lb, ub
   const int len = PyObject_Length(pCutList);
   // loop through each cut
   // We can use these to construct a C++ DecompVar objects
   double lb, ub;
   PyObject* pRow, *pLb, *pUb;
   string name;
   double value;

   for (int i = 0; i < len; i++) {
      pRow = PySequence_GetItem(pCutList, i);
      char getLb[] = "getLb";
      pLb = PyObject_CallMethod(pRow, getLb, NULL);

      if (pLb == NULL) {
         throw UtilException("Error calling method row.getLb()", "generateCuts", "DippyDecompApp");
      }
      char getUb[] = "getUb";
      pUb = PyObject_CallMethod(pRow, getUb, NULL);

      if (pUb == NULL) {
         throw UtilException("Error calling method row.getUb()", "generateCuts", "DippyDecompApp");
      }

      lb = (pLb == Py_None) ? -m_infinity : PyFloat_AsDouble(pLb);
      ub = (pUb == Py_None) ?  m_infinity : PyFloat_AsDouble(pUb);
      int*     varInds = NULL;
      double* varVals = NULL;
      int numPairs = pyColDict_AsPackedArrays(pRow, m_colIndices, &varInds, &varVals);
      assert(numPairs == PyObject_Length(pRow));
      // arrays are now owned by the Cut object
      DippyDecompCut* cut =  new DippyDecompCut(lb, ub, numPairs, varInds, varVals);
      cutList.push_back(cut);
   }

   assert(!PyErr_Occurred());

   return len;
}

/**
 * APPheuristics callback
 *
 * Called by DIP, we interface with Python
 */
int DippyDecompApp::APPheuristics(const double* xhat, const double* origCost, vector<DecompSolution*>& xhatIPFeas)
{
   if (!m_pyHeuristics) {
      return 0;
   }

   PyObject* pSolution = pyTupleList_FromDoubleArray(xhat, m_colList);
   PyObject* pObjective = pyTupleList_FromDoubleArray(origCost, m_colList);
   char solveHeuristics[] = "solveHeuristics";
   char OO[] = "OO";
   PyObject* pSolList = PyObject_CallMethod(m_pProb, solveHeuristics, OO, pSolution, pObjective);

   if (pSolList == NULL) {
      throw UtilException("Error calling method prob.solveHeuristics()", "APPheuristics",
                          "DippyDecompApp");
   }

   // This should never happen, pyHeuristics should be set to false in dippy.py
   if (pSolList == Py_None)
      // method exists, but is not implemented, return 0
   {
      return 0;
   }

   // APPheuristics returns dictionary of (variable, value) pairs
   const int len = PyObject_Length(pSolList);

   // loop through each solution
   for (int i = 0; i < len; i++) {
      pSolution = PyList_GetItem(pSolList, i);
      int*     varInds = NULL;
      double* varVals = NULL;
      int numPairs = pyColDict_AsPackedArrays(pSolution, m_colIndices, &varInds, &varVals);
      assert(numPairs == PyObject_Length(pSolution));
      double* sol = new double[m_numCols];
      UtilFillN(sol, m_numCols, 0.0);

      for (int j = 0; j < numPairs; j++) {
         sol[varInds[j]] = varVals[j];
      }

      xhatIPFeas.push_back(new DecompSolution(m_numCols, sol, origCost));
      delete [] sol;
      delete [] varInds;
      delete [] varVals;
   }

   assert(!PyErr_Occurred());

   return len;
}

/**
 * generateInitVars callback
 *
 * Called by DIP, we interface with Python
 */
int DippyDecompApp::generateInitVars(DecompVarList& initVars)
{
   if (!m_pyInitVars) {
      return 0;
   }

   char generateInitVars[] = "generateInitVars";
   PyObject* pVarList = PyObject_CallMethod(m_pProb, generateInitVars, NULL);

   if (pVarList == NULL) {
      throw UtilException("Error calling method prob.generateInitVars()", "generateInitVars",
                          "DippyDecompApp");
   }

   if (pVarList == Py_None)
      // method exists, but is not implemented, return 0
   {
      return 0;
   }

   int nVars = PyObject_Length(pVarList);
   // generateInitVars returns 2-tuples (index, (cost, dictionary of (variable, value) pairs))
   // We can use these to construct a C++ DecompVar objects
   double cost;
   PyObject* pColDict;

   for (int i = 0; i < nVars; i++) {
      PyObject* pTuple = PyList_GetItem(pVarList, i);
      int whichBlock = m_relaxIndices[PyTuple_GetItem(pTuple, 0)];
      PyObject* pVarTuple = PyTuple_GetItem(pTuple, 1);
      cost   = PyFloat_AsDouble(PyTuple_GetItem(pVarTuple, 0));
      pColDict = PyTuple_GetItem(pVarTuple, 1);
      int*     varInds = NULL;
      double* varVals = NULL;
      DecompVarType varType; 
      
      int numPairs = pyColDict_AsPackedArrays(pColDict, m_colIndices, &varInds, &varVals, varType);
      assert(numPairs == PyObject_Length(pColDict));
      DecompVar* var =  new DecompVar(numPairs, varInds, varVals, cost, varType);
      var->setBlockId(whichBlock);
      initVars.push_back(var);
   }

   assert(!PyErr_Occurred());

   return nVars;
}
