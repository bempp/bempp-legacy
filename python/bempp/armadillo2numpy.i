// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

%include "numpy.i"

%header %{
#include <algorithm> // std::copy
#include <armadillo>
#include <iostream>
#include <new>
%}

// This must be called at the start of each module to import numpy.
%init %{
    import_array();
%}

%define %arma_numpy_typemaps(DATA_TYPE, DATA_TYPECODE)

/************************/
/* Input Array Typemaps */
/************************/

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
    (const arma::Col< DATA_TYPE >& IN_COL)
{
    $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in,
         fragment="NumPy_Fragments")
    (const arma::Col< DATA_TYPE >& IN_COL)
    (PyArrayObject* array=NULL, int is_new_object=0, arma::Col< DATA_TYPE > arma_array)
{
    array = obj_to_array_contiguous_allow_conversion($input, DATA_TYPECODE,
                                                   &is_new_object);
    if (!array || !require_dimensions(array, 1))
        SWIG_fail;
    arma_array = arma::Col< DATA_TYPE >((DATA_TYPE*) array_data(array),
                               array_size(array, 0),
                               false); // don't copy data
    $1 = &arma_array;
}
%typemap(freearg)
    (const arma::Col< DATA_TYPE >& IN_COL)
{
    if (is_new_object$argnum && array$argnum) {
        Py_DECREF(array$argnum); 
    }
}

/* ------------------------------------------------------------------------- */

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
    (const arma::Mat< DATA_TYPE >& IN_MAT)
{
    $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in,
         fragment="NumPy_Fragments")
    (const arma::Mat< DATA_TYPE >& IN_MAT)
    (PyArrayObject* array=NULL, int is_new_object=0, arma::Mat< DATA_TYPE > arma_array)
{
    array = obj_to_array_fortran_allow_conversion($input, DATA_TYPECODE,
        &is_new_object); 
    if (!array || !require_dimensions(array, 2) || !require_fortran(array))
        SWIG_fail;
    arma_array = arma::Mat< DATA_TYPE >((DATA_TYPE*) array_data(array),
                                      array_size(array, 0),
                                      array_size(array, 1),
                                      false); // don't copy data
    $1 = &arma_array;
}
%typemap(freearg)
    (const arma::Mat< DATA_TYPE >& IN_MAT)
{
  if (is_new_object$argnum && array$argnum) {
      Py_DECREF(array$argnum); 
  }
}

/* ------------------------------------------------------------------------- */
/* This typemap behaves like IN_MAT, but it also makes the wrapped function 
return pointers to the SWIG wrappers of the Numpy array and the Armadillo array 
(allocated on the heap). */

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
    (const arma::Mat< DATA_TYPE >& IN_MAT_OUT_WRAPPERS)
{
    $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in, fragment="NumPy_Fragments")
    (const arma::Mat< DATA_TYPE >& IN_MAT_OUT_WRAPPERS)
    (PyArrayObject* array=NULL, int is_new_object=0, arma::Mat< DATA_TYPE >* arma_array)
{
    array = obj_to_array_fortran_allow_conversion($input, DATA_TYPECODE,
        &is_new_object); 
    if (!array || !require_dimensions(array, 2) || !require_fortran(array))
        SWIG_fail;
    arma_array = new arma::Mat< DATA_TYPE >((DATA_TYPE*) array_data(array),
                                          array_size(array, 0),
                                          array_size(array, 1),
                                          false); // don't copy data
    if (!arma_array) {
        if (is_new_object && array)
            Py_DECREF(array);
        SWIG_fail;
    }
    $1 = arma_array;
}
%typemap(argout)
    (const arma::Mat< DATA_TYPE >& IN_MAT_OUT_WRAPPERS)
{
    PyObject* wrapper_arma_mat = 
        SWIG_NewPointerObj(%as_voidptr(arma_array$argnum), 
                           $1_descriptor,
                           SWIG_POINTER_OWN | %newpointer_flags);
    if (!wrapper_arma_mat) {
        if (arma_array$argnum)
            delete arma_array$argnum;
        if (is_new_object$argnum && array$argnum)
            Py_DECREF(array$argnum);
        SWIG_fail;
    }
    $result = SWIG_Python_AppendOutput($result, reinterpret_cast<PyObject*>(array$argnum));
    $result = SWIG_Python_AppendOutput($result, wrapper_arma_mat);
}

/*************************/
/* Output Array Typemaps */
/*************************/
 
%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (arma::Col< DATA_TYPE >& ARGOUT_COL)
    (PyArrayObject* array=NULL, arma::Col< DATA_TYPE > arma_array)
{
    $1 = &arma_array;
}
%typemap(argout)
    (arma::Col< DATA_TYPE >& ARGOUT_COL)
{
    npy_intp dims[1];
    dims[0] = arma_array$argnum.n_rows;
    array$argnum = reinterpret_cast<PyArrayObject*>(
        PyArray_EMPTY(1, dims, DATA_TYPECODE, NPY_FORTRAN));
    if (!array$argnum)
        SWIG_fail;
    std::copy(arma_array$argnum.begin(), arma_array$argnum.end(),  
        reinterpret_cast<DATA_TYPE*>(array_data(array$argnum)));
    $result = SWIG_Python_AppendOutput($result, 
        reinterpret_cast<PyObject*>(array$argnum));
}

/* ------------------------------------------------------------------------- */

%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (arma::Row< DATA_TYPE >& ARGOUT_ROW)
    (PyArrayObject* array=NULL, arma::Row< DATA_TYPE > arma_array)
{
    $1 = &arma_array;
}
%typemap(argout)
    (arma::Row< DATA_TYPE >& ARGOUT_ROW)
{
    npy_intp dims[1];
    dims[0] = arma_array$argnum.n_cols;
    array$argnum = reinterpret_cast<PyArrayObject*>(
        PyArray_EMPTY(1, dims, DATA_TYPECODE, NPY_FORTRAN));
    if (!array$argnum)
        SWIG_fail;
    std::copy(arma_array$argnum.begin(), arma_array$argnum.end(),  
	    reinterpret_cast<DATA_TYPE*>(array_data(array$argnum)));
    $result = SWIG_Python_AppendOutput($result, 
        reinterpret_cast<PyObject*>(array$argnum));
}

/* ------------------------------------------------------------------------- */

%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (arma::Mat< DATA_TYPE >& ARGOUT_MAT)
    (PyArrayObject* array=NULL, arma::Mat< DATA_TYPE > arma_array)
{
    $1 = &arma_array;
}
%typemap(argout)
    (arma::Mat< DATA_TYPE >& ARGOUT_MAT)
{
    npy_intp dims[2];
    dims[0] = arma_array$argnum.n_rows;
    dims[1] = arma_array$argnum.n_cols;
    array$argnum = reinterpret_cast<PyArrayObject*>(
        PyArray_EMPTY(2, dims, DATA_TYPECODE, NPY_FORTRAN));
    if (!array$argnum)
        SWIG_fail;
    std::copy(arma_array$argnum.begin(), arma_array$argnum.end(),  
        reinterpret_cast<DATA_TYPE*>(array_data(array$argnum)));
    $result = SWIG_Python_AppendOutput($result, 
        reinterpret_cast<PyObject*>(array$argnum));
}

/* ------------------------------------------------------------------------- */

%typemap(in, numinputs=0,
         fragment="NumPy_Fragments")
    (arma::Cube< DATA_TYPE >& ARGOUT_CUBE)
    (PyArrayObject* array=NULL, arma::Cube< DATA_TYPE > arma_array)
{
    $1 = &arma_array;
}
%typemap(argout)
    (arma::Cube< DATA_TYPE >& ARGOUT_CUBE)
{
    npy_intp dims[3];
    dims[0] = arma_array$argnum.n_rows;
    dims[1] = arma_array$argnum.n_cols;
    dims[2] = arma_array$argnum.n_slices;
    array$argnum = reinterpret_cast<PyArrayObject*>(
        PyArray_EMPTY(3, dims, DATA_TYPECODE, NPY_FORTRAN));
    if (!array$argnum)
        SWIG_fail;
    std::copy(arma_array$argnum.begin(), arma_array$argnum.end(),  
        reinterpret_cast<DATA_TYPE*>(array_data(array$argnum)));
    $result = SWIG_Python_AppendOutput($result, 
        reinterpret_cast<PyObject*>(array$argnum));
}

/**************************/
/* Inplace Array Typemaps */
/**************************/

/* Typemap suite for (DATA_TYPE* INPLACE_ARRAY1, DIM_TYPE DIM1)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (arma::Col< DATA_TYPE >& INPLACE_COL)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),
                                                 DATA_TYPECODE);
}
%typemap(in,
         fragment="NumPy_Fragments")
  (arma::Col< DATA_TYPE >& INPLACE_COL)
  (PyArrayObject* array=NULL, arma::Col< DATA_TYPE > arma_array)
{
  array = obj_to_array_no_conversion($input, DATA_TYPECODE);
  if (!array || !require_dimensions(array,1) || !require_contiguous(array)
      || !require_native(array)) SWIG_fail;
  // Use placement new to reinitialise the Armadillo array using the
  // "advanced" constructor taking a pointer to existing data. This is needed
  // because SWIG initialises variables with the default constructor.
  // (Another way would be to allocate a new Col object on the heap),
  arma_array.~Col< DATA_TYPE >();
  new (&arma_array) arma::Col< DATA_TYPE >((DATA_TYPE*) array_data(array),
                                           array_size(array, 0),
                                           false); // don't copy data
  $1 = &arma_array;
}

%enddef    /* %arma_numpy_typemaps() macro */
/* *************************************************************** */

/* Concrete instances of the %arma_numpy_typemaps() macro: Each invocation
 * below applies all of the typemaps above to the specified data type.
 */
%arma_numpy_typemaps(signed char       , NPY_BYTE     )
%arma_numpy_typemaps(char              , NPY_CHAR     )
%arma_numpy_typemaps(unsigned char     , NPY_UBYTE    )
%arma_numpy_typemaps(short             , NPY_SHORT    )
%arma_numpy_typemaps(unsigned short    , NPY_USHORT   )
%arma_numpy_typemaps(int               , NPY_INT      )
%arma_numpy_typemaps(unsigned int      , NPY_UINT     )
%arma_numpy_typemaps(long              , NPY_LONG     )
%arma_numpy_typemaps(unsigned long     , NPY_ULONG    )
%arma_numpy_typemaps(long long         , NPY_LONGLONG )
%arma_numpy_typemaps(unsigned long long, NPY_ULONGLONG)
%arma_numpy_typemaps(float             , NPY_FLOAT    )
%arma_numpy_typemaps(double            , NPY_DOUBLE   )

%arma_numpy_typemaps(Bempp::ctype      , NPY_DOUBLE   )

%arma_numpy_typemaps(std::complex<float>,  NPY_CFLOAT  );
%arma_numpy_typemaps(std::complex<double>, NPY_CDOUBLE );
