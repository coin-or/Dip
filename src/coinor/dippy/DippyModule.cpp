#include "Python.h"
#include "Decomp.h"

#ifdef _WIN32
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT extern "C"
#endif

// prototype
DLLEXPORT PyObject* Solve(PyObject* self, PyObject* args);

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

// methods exposed by this module
static PyMethodDef dippy_module_methods[] = {
   {"Solve", Solve, METH_VARARGS, "Solve"},

   {NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3

static int dippy_module_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int dippy_module_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_dippy",
        NULL,
        sizeof(struct module_state),
        dippy_module_methods,
        NULL,
        dippy_module_traverse,
        dippy_module_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit__dippy(void)

#else
#define INITERROR return

void
init_dippy(void)
#endif
{
  PyObject *module, *d, *s;
#if PY_MAJOR_VERSION >= 3
  module=PyModule_Create(&moduledef);
#else
  module=Py_InitModule("_dippy", dippy_module_methods);
#endif
  d = PyModule_GetDict(module);
  s = PyUnicode_FromString("0.2");
  PyDict_SetItemString(d, "__version__", s);
  s = PyUnicode_FromString("See polyhedron.py");
  PyDict_SetItemString(d, "__doc__", s);
  //dippy_error = PyUnicode_FromString("_dippy.error");
  //PyDict_SetItemString(d, "error", cdd_error);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _cdd");
#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}


