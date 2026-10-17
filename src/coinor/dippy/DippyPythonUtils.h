#ifndef DIPPY_PYTHONUTILS_INCLUDED
#define DIPPY_PYTHONUTILS_INCLUDED

#include "Python.h"

#include "Decomp.h"
#include "DecompAlgo.h"

#include <map>
#include <vector>
using namespace std;

// Some convenience functions for converting between Python objects and
// C/C++ data structures

/**
 * Convert a double array to a Python tuple list
 *
 * The list is (Python object, value) tuples
 *
 * Returns Python tuple list with length = pList length
 */
PyObject* pyTupleList_FromDoubleArray(const double* values, PyObject* pList);

/**
 * Package a AlpsDecompTreeNode using a DecompAlgo into a
 * Python list
 *
 * The list is (Python object, value) tuples
 *
 * Returns Python tuple list
 */
PyObject* pyTupleList_FromNode(DecompAlgo* algo, DecompStatus decompStatus);

/**
 * Convert a column dictionary to a (int, double) vector
 *
 * The dictionary has columns as keys
 * and coefficients as values
 *
 */
void pyColDict_AsPairedVector(PyObject* pColDict, vector< pair<int, double> >& vector, map<PyObject*, int> indices);

/**
 * Convert a column dictionary to packed arrays
 *
 * The dictionary has columns as keys
 * and coefficients as values
 *
 * Returns length of index and value arrays
 */
int pyColDict_AsPackedArrays(PyObject* pColDict, map<PyObject*, int> indices, int** inds, double** vals);

int pyColDict_AsPackedArrays(PyObject* pColDict, map<PyObject*, int> indices, int** inds, double** vals, DecompVarType & varType);

/**
 * Convert a list of Python constraints to a CoinPackedMatrix
 *
 * The format of each constraint is as a dictionary with variables as keys
 * and coefficients as values
 */

CoinPackedMatrix* pyConstraints_AsPackedMatrix(PyObject* pRowList,
      map<PyObject*, int> rowIndices, map<PyObject*, int> colIndices);

/**
 * Creates a (key,value) tuple and appends to a Python list of tuples *
 */
void addTupleToPyList(PyObject* pList, PyObject* key, PyObject* value);

/**
 * Creates a (key,value) tuple and inserts in a Python list of tuples *
 */
void insertTupleToPyList(PyObject* pList, unsigned position, PyObject* key, PyObject* value);

#endif
