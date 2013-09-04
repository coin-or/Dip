#include "Python.h"
#include "Decomp.h"

#ifdef _WIN32
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT extern "C"
#endif

// prototype
DLLEXPORT PyObject* Solve(PyObject* self, PyObject* args);

// methods exposed by this module
static PyMethodDef Methods[] = {
   {"Solve", Solve, METH_VARARGS, "Solve"},

   {NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC
init_dippy(void)
{
   PyObject* pMod = Py_InitModule("_dippy", Methods);
}


