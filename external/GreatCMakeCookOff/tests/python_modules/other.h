#ifndef TEST_PYTHON_MODULE
#define TEST_PYTHON_MODULE

#include <Python.h>

#ifdef __cplusplus
    extern "C" 
#endif
PyObject * meaning_of_life(PyObject *_module);

#endif
