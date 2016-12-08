#include <Python.h>
#include "other.h"

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#   define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#   define GETSTATE(m) (&_state)
    static struct module_state _state;
#endif

static PyObject * error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef extension_methods[] = {
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {"meaning_of_life", (PyCFunction)meaning_of_life, METH_NOARGS, NULL},
    {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

    static int extension_traverse(PyObject *m, visitproc visit, void *arg) {
        Py_VISIT(GETSTATE(m)->error);
        return 0;
    }

    static int extension_clear(PyObject *m) {
        Py_CLEAR(GETSTATE(m)->error);
        return 0;
    }


    static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "extension",
            NULL,
            sizeof(struct module_state),
            extension_methods,
            NULL,
            extension_traverse,
            extension_clear,
            NULL
    };

#   define INITERROR return NULL

    PyObject * PyInit_extension(void)

#else
#   define INITERROR return

    void initextension(void)
#endif

{
#   if PY_MAJOR_VERSION >= 3
        PyObject *module = PyModule_Create(&moduledef);
#   else
        PyObject *module = Py_InitModule("extension", extension_methods);
#   endif

    if (module == NULL) INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("extension.Error", PyExc_Exception, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    PyModule_AddObject(module, "Error", st->error);

#   if PY_MAJOR_VERSION >= 3
        return module;
#   endif
}
