#include "DippyPythonUtils.h"
#include "DippyDecompApp.h"
#include "UtilMacros.h"
#include "UtilMacrosDecomp.h"
#include "AlpsDecompNodeDesc.h"

#include <assert.h>

//===========================================================================//

// Some convenience functions for converting between Python objects and
// C/C++ data structures

/**
 * Convert a double array to a Python tuple list
 *
 * The list is (Python object, value) tuples
 *
 * Returns Python tuple list with length = pList length
 */
PyObject* pyTupleList_FromDoubleArray(const double* values, PyObject* pList)
{
   int len = PyObject_Length(pList);
   PyObject* pTupleList = PyList_New(len), *pObj;

   for (int i = 0; i < len; i++) {
      pObj = PyList_GetItem(pList, i);
      Py_XINCREF(pObj);
      insertTupleToPyList(pTupleList, i, pObj, PyFloat_FromDouble(values[i]));
   }

   assert(!PyErr_Occurred());

   return pTupleList;
}

/**
 * Package a AlpsDecompTreeNode using a DecompAlgo into a
 * Python list
 *
 * The list is (Python object, value) tuples
 *
 * Returns Python tuple list
 */
PyObject* pyTupleList_FromNode(DecompAlgo* algo, DecompStatus decompStatus)
{
   PyObject* pOutput = PyList_New(0);
   AlpsDecompTreeNode* node = (AlpsDecompTreeNode*) algo->getCurrentNode();
   double lb = algo->getObjBestBoundLB(), ub = algo->getObjBestBoundUB();
   double quality = node->getQuality();
   string status;

   switch (decompStatus) {
   case STAT_IP_FEASIBLE:
      status = "Solution";
      break;

   case STAT_FEASIBLE:
      if (lb > quality) {
         quality = lb;
      }

      if (quality >= ub) {
         status = "Pruned";
      } else {
         status = "Candidate";
      }

      break;

   case STAT_INFEASIBLE:
      status = "Infeasible";
      break;

   default:
      status = "Unknown";
   }

   // Add into the list the output needed
   addTupleToPyList(pOutput, PyUnicode_FromString("nodeIndex"),
                    PyLong_FromLong(node->getIndex()));
   addTupleToPyList(pOutput, PyUnicode_FromString("parentIndex"),
                    PyLong_FromLong(node->getParentIndex()));
   addTupleToPyList(pOutput, PyUnicode_FromString("nodeDepth"),
                    PyLong_FromLong(node->getDepth()));
   addTupleToPyList(pOutput, PyUnicode_FromString("nodeQuality"),
                    PyFloat_FromDouble(quality));
   addTupleToPyList(pOutput, PyUnicode_FromString("globalLB"),
                    PyFloat_FromDouble(lb));
   addTupleToPyList(pOutput, PyUnicode_FromString("globalUB"),
                    PyFloat_FromDouble(ub));
   addTupleToPyList(pOutput, PyUnicode_FromString("nodeStatus"),
                    PyUnicode_FromString(status.c_str()));
   addTupleToPyList(pOutput, PyUnicode_FromString("branchedDir"),
                    PyLong_FromLong(dynamic_cast<AlpsDecompNodeDesc*>
                                   (algo->getCurrentNode()->getDesc())->getBranchedDir()));
   // Copy the current solution into a Python list
   const double* xhat = algo->getXhat();
   DippyDecompApp* app = (DippyDecompApp*)algo->getDecompApp();
   PyObject* pSolutionList = pyTupleList_FromDoubleArray(xhat,
                             app->m_colList);
   addTupleToPyList(pOutput, PyUnicode_FromString("xhat"), pSolutionList);
   /** MO. 29/2/2012 - This section was originally an attempt to add "simple" cuts. i.e.,
       that the sum of non-basic variables >= 1 (or at least a variant for lb and ub), so I
   	passed the variable bounds at the node and lists of basic, lower bound and upper bound
   	variables from the original problem. However, the cuts added "clashed" with the CGL cuts,
   	probably because the variables introduced by these cuts, e.g., slacks, were not considered.
   	I have abandoned this direction for now, but I like the idea of getting the full node problem
   	so more complex cuts that use, e.g., the basis, can be implemented in Python.

   int numOrig = algo->getModelCore().getModel()->getNumCols();

   const double * lb = algo->getColLBNode();
   const double * ub = algo->getColUBNode();
   PyObject * pBoundList = PyList_New(0),
   	     * pExtraBoundList = PyList_New(0);
   PyObject * pBoundPair;

   for (int j=0; j<numOrig; j++) {
   	pBoundPair = PyTuple_New(2);
   	if (lb[j] <= -DecompInf)
   		PyTuple_SetItem(pBoundPair, 0, Py_None);
   	else
   		PyTuple_SetItem(pBoundPair, 0, PyFloat_FromDouble(lb[j]));
   	if (ub[j] >= DecompInf)
   		PyTuple_SetItem(pBoundPair, 1, Py_None);
   	else
   		PyTuple_SetItem(pBoundPair, 1, PyFloat_FromDouble(ub[j]));
   	addTupleToPyList(pBoundList, PyList_GetItem(app->m_colList, j), pBoundPair);
   }
   addTupleToPyList(pOutput, PyUnicode_FromString("bounds"), pBoundList);

   // Copy the original variables into "status" lists
   PyObject * pBasisList = PyList_New(0),
   	     * pLBList    = PyList_New(0),
   		 * pUBList    = PyList_New(0);

   if (algo->getMasterOSI()->basisIsAvailable()) {
   	int numRows = algo->getMasterOSI()->getNumRows(),
   		numCols = algo->getMasterOSI()->getNumCols();
   	int * rstat = new int[numRows],
   		* cstat = new int [numCols];

   	algo->getMasterOSI()->getBasisStatus(cstat, rstat);
   	// MO (28/2/2012) - Assuming that any extra variables are added at the end, is this true?
   	for (int j=0; j<numOrig; j++)
   		switch (cstat[j]) {
   		case 1: // basic
   			PyList_Append(pBasisList, PyList_GetItem(app->m_colList, j));
   			break;
   		case 2: // upper
   			PyList_Append(pUBList, PyList_GetItem(app->m_colList, j));
   			break;
   		case 3: // lower
   			PyList_Append(pLBList, PyList_GetItem(app->m_colList, j));
   			break;
   		default:
             throw UtilException("Error calling method pyTupleList_FromNode()", "pyTupleList_FromNode", "DippyPythonUtils");
   		}

   	delete [] rstat;
   	delete [] cstat;
   }

   addTupleToPyList(pOutput, PyUnicode_FromString("basic"), pBasisList);
   addTupleToPyList(pOutput, PyUnicode_FromString("atLB"), pLBList);
   addTupleToPyList(pOutput, PyUnicode_FromString("atUB"), pUBList);
   */

   assert(!PyErr_Occurred());

   return pOutput;
}

/**
 * Convert a column dictionary to a (int, double) vector
 *
 * The dictionary has columns as keys
 * and coefficients as values
 *
 */
void pyColDict_AsPairedVector(PyObject* pColDict, vector<pair<int, double> >& vec, map<PyObject*, int> indices)
{
   int len = PyObject_Length(pColDict);
   vec.clear();
   PyObject* pKeys = PyDict_Keys(pColDict), *pCol;
   double value;
   int index;

   for (int i = 0; i < len; i++) {
      pCol = PyList_GetItem(pKeys, i);
      value = PyFloat_AsDouble(PyDict_GetItem(pColDict, pCol));
      index = indices[pCol];

      if ( (index < 0) || (index >= indices.size()) ) {
	 char str[] = "__str__";
	 PyObject* pColName = PyObject_CallMethod(pCol, str, NULL);

         if (pColName == NULL) {
            throw UtilException("Error calling method col.__str__()", "pyColDict_AsPairedVector",
                                "DippyPythonUtils");
         }

         string name = PyBytes_AsString(PyUnicode_AsEncodedString(pColName, "UTF-8", "strict"));

         throw UtilException("Bad index for " + name, "pyTupleList_AsPairedVector",
                             "DippyPythonUtils");
      }

      vec.push_back(pair<int, double>(index, value));
   }

   assert(!PyErr_Occurred());
}

/**
 * Convert a column dictionary to packed arrays
 *
 * The dictionary has columns as keys
 * and coefficients as values
 *
 * Returns length of index and value arrays
 */
int pyColDict_AsPackedArrays(PyObject* pColDict, map<PyObject*, int> indices, int** inds, double** vals)
{
   int len = PyObject_Length(pColDict);
   *inds = new int[len];
   *vals = new double[len];
   PyObject* pKeys = PyDict_Keys(pColDict);
   PyObject* pCol;
   double value;
   int index;

   for (int i = 0; i < len; i++) {
      pCol = PyList_GetItem(pKeys, i);
      value = PyFloat_AsDouble(PyDict_GetItem(pColDict, pCol));
      index = indices[pCol];

      if ( (index < 0) || (index >= indices.size()) ) {
	char getName[] = "getName";
         PyObject* pColName = PyObject_CallMethod(pCol, getName, NULL);

         if (pColName == NULL) {
            throw UtilException("Error calling method col.getName()", "pyColDict_AsPackedArrays",
                                "DippyPythonUtils");
         }

         string name = PyBytes_AsString(PyUnicode_AsEncodedString(pColName, "UTF-8", "strict"));

         throw UtilException("Bad index for " + name, "pyColDict_AsPackedArrays", "DippyPythonUtils");
      }

      (*inds)[i] = index;
      (*vals)[i] = value;
   }

   assert(!PyErr_Occurred());

   return len;
}

int pyColDict_AsPackedArrays(PyObject* pColDict, map<PyObject*, int> indices, int** inds, double** vals, DecompVarType& varType)
{
   int len = PyObject_Length(pColDict);
   *inds = new int[len];
   *vals = new double[len];
   PyObject* pKeys = PyDict_Keys(pColDict);
   PyObject* pCol;
   double value;
   int index;

   for (int i = 0; i < len; i++) {
      pCol = PyList_GetItem(pKeys, i);
      value = PyFloat_AsDouble(PyDict_GetItem(pColDict, pCol));
      index = indices[pCol];

      if ( (index < 0) || (index >= indices.size()) ) {
	char getName[] = "getName";
         PyObject* pColName = PyObject_CallMethod(pCol, getName, NULL);

         if (pColName == NULL) {
            throw UtilException("Error calling method col.getName()", "pyColDict_AsPackedArrays",
                                "DippyPythonUtils");
         }

         string name = PyBytes_AsString(PyUnicode_AsEncodedString(pColName, "UTF-8", "strict"));

         throw UtilException("Bad index for " + name, "pyColDict_AsPackedArrays",
                             "DippyPythonUtils");
      }
      char getVarType[] = "getVarType";
      PyObject* pColType = PyObject_CallMethod(pCol, getVarType, NULL);      
      if (pColType == NULL){
         throw UtilException("getVarType call failed.", "pyColDict_AsPackedArrays",
                             "DippyPythonUtils");
      }
     
      (*inds)[i] = index;
      (*vals)[i] = value;
   }

   assert(!PyErr_Occurred());

   return len;
}

/**
 * Convert a list of Python constraints to a CoinPackedMatrix
 *
 * The format of each constraint is as a dictionary with variables as keys
 * and coefficients as values
 */

CoinPackedMatrix* pyConstraints_AsPackedMatrix(PyObject* pRowList,
      map<PyObject*, int> rowIndices, map<PyObject*, int> colIndices)
{
   int len = PyObject_Length(pRowList);
   int rowInd, colInd, num;
   string rowName, colName;
   double val, lb, ub;
   PyObject* pRow, *pKeys, *pCol;
   // First get the total number of non-zeros from all the tuples
   int numNZs = 0;

   for (int i = 0; i < len; i++) {
      pRow = PyList_GetItem(pRowList, i);
      num = PyObject_Length(pRow);
      numNZs += num;
   }

   // Now read and process the tuples
   int start = 0;
   int* rowInds = new int[numNZs];
   UtilFillN(rowInds, numNZs, -1);
   int* colInds = new int[numNZs];
   UtilFillN(colInds, numNZs, -1);
   double* values = new double[numNZs];
   UtilFillN(values, numNZs, 0.0);

   for (int i = 0; i < len; i++) {
      pRow = PyList_GetItem(pRowList, i);
      rowInd = rowIndices[pRow];
      pKeys = PyDict_Keys(pRow);
      num = PyObject_Length(pKeys);

      for (int n = 0; n < num; n++) {
         pCol = PyList_GetItem(pKeys, n);
         colInd = colIndices[pCol];
         val = PyFloat_AsDouble(PyDict_GetItem(pRow, pCol));
         rowInds[start + n] = rowInd;
         colInds[start + n] = colInd;
         values[start + n] = val;
      }

      start += num;
   }

   assert(!PyErr_Occurred());

   return new CoinPackedMatrix(false, rowInds, colInds, values, numNZs);
}

/**
 * Creates a (key,value) tuple and appends to a Python list of tuples *
 */
void addTupleToPyList(PyObject* pList, PyObject* key, PyObject* value)
{
   PyObject* pTuple = PyTuple_New(2);
   PyTuple_SetItem(pTuple, 0, key);
   PyTuple_SetItem(pTuple, 1, value);
   PyList_Append(pList, pTuple);
   assert(!PyErr_Occurred());
}

/**
 * Creates a (key,value) tuple and inserts in a Python list of tuples *
 */
void insertTupleToPyList(PyObject* pList, unsigned position, PyObject* key, PyObject* value)
{
   PyObject* pTuple = PyTuple_New(2);
   PyTuple_SetItem(pTuple, 0, key);
   PyTuple_SetItem(pTuple, 1, value);
   PyList_SetItem(pList, position, pTuple);
   assert(!PyErr_Occurred());
}

