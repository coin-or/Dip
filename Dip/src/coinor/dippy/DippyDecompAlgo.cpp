#include "Python.h"

#include "DippyDecompAlgo.h"
#include "DippyDecompApp.h"
#include "DippyPythonUtils.h"

/**
 * Perform branching
 *
 * This function should populate the (down|up)Branch(LB|UB) vectors with (indices, bound-value) pairs
 * which define the branch.
 */
bool DippyAlgoMixin::chooseBranchSet(DecompAlgo* algo,
				     std::vector< std::pair<int, double> >& downBranchLB,
                                     std::vector< std::pair<int, double> >& downBranchUB,
                                     std::vector< std::pair<int, double> >& upBranchLB,
                                     std::vector< std::pair<int, double> >& upBranchUB)
{
   bool ret_val;

   if (!m_utilParam->GetSetting("pyBranchMethod", true)) {
      return algo->DecompAlgo::chooseBranchSet(downBranchLB, downBranchUB, 
					       upBranchLB, upBranchUB);
   }

   DippyDecompApp* app = (DippyDecompApp*)algo->getDecompApp();
   // copy the current solution into a Python list
   const double* xhat = algo->getXhat();
   PyObject* pSolutionList = pyTupleList_FromDoubleArray(xhat, app->m_colList);
   // try to call chooseBranchSet on the DipProblem python object
   char arg1[] = "chooseBranchSet";
   char arg2[] = "O";
   PyObject* pResult = PyObject_CallMethod(m_pProb, arg1, arg2,  pSolutionList);

   if (pResult == NULL) {
      //something's gone wrong with the function call, a Python exception has 
      //been set 
      throw UtilException("Error calling method prob.chooseBranchSet()", 
			  "chooseBranchSet", "DippyDecompAlgo");
   }

   // need more error checking/assertion setting here
   if (pResult == Py_None) {
      // if chooseBranchSet returns None, do default branching for this algo
      ret_val = algo->DecompAlgo::chooseBranchSet(downBranchLB, downBranchUB, upBranchLB, upBranchUB);
      
      // Original comment: No branching set was returned. This shouldn't happen
      // tkr 11/11/15: Actually, it can happen if the solution is integral, but not feasible.
      // This happens sometimes when column generation is halted because of tailoff and
      // the solution to the relaxation is feasible. I'm leaving the commented code here for
      // posterity
      //assert(ret_val == true);

      //if (!ret_val){
      // throw UtilException("No branch set found in prob.chooseBranchSet()", 
      //		     "chooseBranchSet", "DippyDecompAlgo");
      //}

      if (downBranchUB.size() > 0) {
         PyObject* downBranchVar, * upBranchVar;
         pDownLB = PyDict_New(); // Down branch LBs is an empty dictionary
         pDownUB = PyDict_New();
         downBranchVar = PyList_GetItem(app->m_colList, downBranchUB[0].first);
         PyDict_SetItem(pDownUB, downBranchVar, 
			PyLong_FromLong(static_cast<int>(round(downBranchUB[0].second))));
         pUpLB = PyDict_New();
         upBranchVar = PyList_GetItem(app->m_colList, upBranchLB[0].first);
         PyDict_SetItem(pUpLB, upBranchVar, 
			PyLong_FromLong(static_cast<int>(round(upBranchLB[0].second))));
         pUpUB = PyDict_New(); // Up branch UBs is an empty dictionary
         assert(downBranchVar == upBranchVar);
      }else{
	 //No branching set was returned. Zero out pointers to old branching 
	 //sets
	 assert(ret_val == false);
	 pDownLB = NULL;
	 pDownUB = NULL;
	 pUpLB = NULL;
	 pUpUB = NULL;
      }
      return ret_val;
   } else {
      // or else, the function should have returned 4 lists of (name, value) tuples
      pDownLB = PyTuple_GetItem(pResult, 0);
      pDownUB = PyTuple_GetItem(pResult, 1);
      pUpLB = PyTuple_GetItem(pResult, 2);
      pUpUB = PyTuple_GetItem(pResult, 3);
      // copy the python dictionaries into the result vectors
      pyColDict_AsPairedVector(pDownLB, downBranchLB, app->m_colIndices);
      pyColDict_AsPairedVector(pDownUB, downBranchUB, app->m_colIndices);
      pyColDict_AsPairedVector(pUpLB, upBranchLB, app->m_colIndices);
      pyColDict_AsPairedVector(pUpUB, upBranchUB, app->m_colIndices);
      return true;
   }
   assert(!PyErr_Occurred());
}

void DippyAlgoMixin::postProcessBranch(DecompAlgo* algo, 
				       DecompStatus decompStatus)
{
   PyObject* pOutput = PyList_New(0);

   if (!m_utilParam->GetSetting("pyPostProcessBranch", true)) {
      return;
   }

   AlpsDecompTreeNode* node = (AlpsDecompTreeNode*)algo->getCurrentNode();
   double quality = node->getQuality();

   if (pDownLB != NULL) {
      addTupleToPyList(pOutput, PyUnicode_FromString("pDownLB"), pDownLB);
   }

   if (pDownUB != NULL) {
      addTupleToPyList(pOutput, PyUnicode_FromString("pDownUB"), pDownUB);
   }

   if (pUpLB != NULL) {
      addTupleToPyList(pOutput, PyUnicode_FromString("pUpLB"), pUpLB);
   }

   if (pUpUB != NULL) {
      addTupleToPyList(pOutput, PyUnicode_FromString("pUpUB"), pUpUB);
   }

   addTupleToPyList(pOutput, PyUnicode_FromString("nodeIndex"), PyLong_FromLong(algo->getNodeIndex()));
   addTupleToPyList(pOutput, PyUnicode_FromString("nodeQuality"), PyFloat_FromDouble(quality));
   char arg1[] = "postProcessBranch";
   char arg2[] = "O";
   PyObject* pResult = PyObject_CallMethod(m_pProb, arg1, arg2, pOutput);
   if (pResult == NULL){
      throw UtilException("postProcessNode call failed.", "postProcessNode", "DippyAlgoMixin");
   }

   assert(!PyErr_Occurred());
}

void DippyAlgoMixin::postProcessNode(DecompAlgo* algo, DecompStatus decompStatus)
{
   if (!m_utilParam->GetSetting("pyPostProcessNode", true)) {
      return;
   }

   PyObject* pOutput = pyTupleList_FromNode(algo, decompStatus);
   char arg1[] = "postProcessNode";
   char arg2[] = "O";
   PyObject* pResult = PyObject_CallMethod(m_pProb, arg1, arg2, pOutput);
   if (pResult == NULL){
      throw UtilException("postProcessNode call failed.", "postProcessNode", "DippyAlgoMixin");
   }

   assert(!PyErr_Occurred());
}

