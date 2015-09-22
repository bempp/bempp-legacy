#ifndef BEMPP_PYTHON_DISCRETE_OPERATOR_CONVERSION_HPP
#define BEMPP_PYTHON_DISCRETE_OPERATOR_CONVERSION_HPP

#include "bempp/core/utils/py_types.hpp"
#include "bempp/assembly/discrete_sparse_boundary_operator.hpp"
#include "bempp/assembly/discrete_dense_boundary_operator.hpp"
#include "bempp/assembly/discrete_hmat_boundary_operator.hpp"
#include "bempp/hmat/hmatrix.hpp"
#include "fiber/scalar_traits.hpp"
#include <Python.h>
#include "numpy/arrayobject.h"
#include <boost/variant.hpp>
#include <type_traits>

namespace Bempp {

template <typename ValueType>
PyObject *py_get_sparse_from_discrete_operator(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  shared_ptr<const RealSparseMatrix> sparseOperator =
      (static_pointer_cast<const DiscreteSparseBoundaryOperator<ValueType>>(
           op.get()))
          ->sparseMatrix();

  PyObject *data;
  PyObject *rowind;
  PyObject *col_ptr;
  PyObject *pModule;
  PyObject *pDict;
  PyObject *pFunc;
  PyObject *pArgsCsc;
  PyObject *pArgsShape;
  PyObject *pArgs;
  PyObject *pM;
  PyObject *pN;
  PyObject *pCsc;

  int M;
  int N;

  // Get the CSC matrix function

  pModule = PyImport_ImportModule("scipy.sparse");
  pDict = PyModule_GetDict(pModule);
  pFunc = PyDict_GetItemString(pDict, "csc_matrix");

  // Now get the parameter objects

  RealSparseMatrix::Index *indexOffset =
      const_cast<RealSparseMatrix::Index *>(sparseOperator->outerIndexPtr());
  RealSparseMatrix::Index *indices =
      const_cast<RealSparseMatrix::Index *>(sparseOperator->innerIndexPtr());
  double *values = const_cast<double *>(sparseOperator->valuePtr());

  npy_intp num_nonzeros = sparseOperator->nonZeros();
  npy_intp num_col_ptr = 1 + sparseOperator->cols();

  M = sparseOperator->rows();
  N = sparseOperator->cols();

  pM = PyInt_FromLong(M);
  pN = PyInt_FromLong(N);

  pArgsShape = PyTuple_New(2);
  PyTuple_SetItem(pArgsShape, 0, pM);
  PyTuple_SetItem(pArgsShape, 1, pN);

  data = PyArray_SimpleNew(1, &num_nonzeros, NPY_FLOAT64);
  rowind = PyArray_SimpleNew(1, &num_nonzeros, NPY_INT);
  col_ptr = PyArray_SimpleNew(1, &num_col_ptr, NPY_INT);

  for (npy_intp i = 0; i < num_nonzeros; ++i)
    *((double *)PyArray_GETPTR1(data, i)) = values[i];

  for (npy_intp i = 0; i < num_nonzeros; ++i)
    *((int *)PyArray_GETPTR1(rowind, i)) = indices[i];

  for (npy_intp i = 0; i < num_col_ptr; ++i)
    *((int *)PyArray_GETPTR1(col_ptr, i)) = indexOffset[i];

  pArgsCsc = PyTuple_New(3);
  pArgs = PyTuple_New(2);

  PyTuple_SetItem(pArgsCsc, 0, data);
  PyTuple_SetItem(pArgsCsc, 1, rowind);
  PyTuple_SetItem(pArgsCsc, 2, col_ptr);

  PyTuple_SetItem(pArgs, 0, pArgsCsc);
  PyTuple_SetItem(pArgs, 1, pArgsShape);

  // Now create the CSC matrix object

  pCsc = PyObject_CallObject(pFunc, pArgs);

  // Now clean up

  Py_DECREF(pModule);
  Py_DECREF(pArgs);

  return pCsc;
}

template <typename ValueType>
PyObject *py_array_from_dense_operator(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {

  return (*static_pointer_cast<const DiscreteDenseBoundaryOperator<ValueType>,
                               const DiscreteBoundaryOperator<ValueType>>(op))
      .asNumpyObject();
}

template <typename ValueType>
shared_ptr<const hmat::HMatrix<ValueType, 2>> py_hmatrix_from_discrete_operator(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {

  return static_pointer_cast<const DiscreteHMatBoundaryOperator<ValueType>>(op)
      ->hMatrix();
}
}

#endif
