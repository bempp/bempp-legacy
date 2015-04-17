// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "bempp/common/config_ahmed.hpp"

#include <iostream>

#include "discrete_boundary_operator.hpp"

#include "discrete_boundary_operator_sum.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "scaled_discrete_boundary_operator.hpp"
#include "transposed_discrete_boundary_operator.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Bempp {

template <typename ValueType>
Matrix<ValueType> DiscreteBoundaryOperator<ValueType>::asMatrix() const {
  // Default brute-force implementation: apply operator to all basis vectors
  const size_t nRows = rowCount();
  const size_t nCols = columnCount();
  Vector<ValueType> unit(nCols);
  Matrix<ValueType> result(nRows, nCols);
  result.setZero(); // for safety, in case there was a bug in the handling of
                   // beta == 0. in a particular subclass' applyBuiltInImpl()
                   // override...
  unit.setZero();
  for (size_t i = 0; i < nCols; ++i) {
    Vector<ValueType> activeCol(result.rows());
    // arma::Col<ValueType> activeCol(result.unsafe_col(i));
    if (i > 0)
      unit(i - 1) = 0.;
    unit(i) = 1.;
    applyBuiltInImpl(NO_TRANSPOSE, unit, activeCol, 1., 0.);
    result.col(i) = activeCol;
  }

  return result;
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::apply(
    const TranspositionMode trans, const Matrix<ValueType> &x_in,
    Matrix<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  bool transposed = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);
  if (x_in.rows() != (transposed ? rowCount() : columnCount()))
    throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                "vector x_in has invalid length");
  if (y_inout.rows() != (transposed ? columnCount() : rowCount()))
    throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                "vector y_inout has invalid length");
  if (x_in.cols() != y_inout.cols())
    throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                "vectors x_in and y_inout must have "
                                "the same number of columns");

  for (size_t i = 0; i < x_in.cols(); ++i) {

      applyBuiltInImpl(trans,Eigen::Ref<Vector<ValueType>>(const_cast<Matrix<ValueType>&>(x_in).col(i)),
                       Eigen::Ref<Vector<ValueType>>(const_cast<Matrix<ValueType>&>(y_inout).col(i)),
                       alpha,beta);
  }
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::apply(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {

    this->apply(trans,
                Eigen::Ref<Vector<ValueType>>(const_cast<Vector<ValueType>&>(x_in)),
                Eigen::Ref<Vector<ValueType>>(const_cast<Vector<ValueType>&>(y_inout)),
                alpha,beta);

}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::apply(const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
           Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
           const ValueType beta) const {

    applyBuiltInImpl(trans, x_in, y_inout, alpha, beta);

}

template <typename ValueType>
PyObject* DiscreteBoundaryOperator<ValueType>::apply(const TranspositionMode trans, const PyObject* x_in) const {

  PyObject* x = const_cast<PyObject*>(x_in);

  Py_INCREF(x); // Increase reference count while in function

  if (!PyArray_Check(x)){
      Py_DECREF(x);
      throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
              "Python object is not of PyArray type");
  }

  int typenum = PyArray_TYPE(reinterpret_cast<PyArrayObject*>(x));
  if (typenum!=Fiber::ScalarTraits<ValueType>::NumpyTypeNum){
      Py_DECREF(x);
      throw std::invalid_argument("DiscreteBoundaryOperator::apply() "
              "Python Array has wrong type.");

  }

  if (!PyArray_ISALIGNED(reinterpret_cast<PyArrayObject*>(x))){
      Py_DECREF(x);
      throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
              "Python array must be aligned.");
  }

  PyObject* x_f;

  if (!PyArray_IS_F_CONTIGUOUS(reinterpret_cast<PyArrayObject*>(x)))
      x_f = PyArray_NewCopy(reinterpret_cast<PyArrayObject*>(x),NPY_FORTRANORDER);
  else
      x_f = x;

  int ndim = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(x_f));

  if (ndim != 2){
      Py_DECREF(x_f);
      throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
              "PyArray x_in must have two dimensions.");
  }

  npy_intp nrows = PyArray_DIM(reinterpret_cast<PyArrayObject*>(x_f),0);
  npy_intp ncols = PyArray_DIM(reinterpret_cast<PyArrayObject*>(x_f),1);

  bool transposed = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);
  if (nrows != (transposed ? rowCount() : columnCount())){
    Py_DECREF(x_f);
    throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                "PyArray x_in has invalid number of rows");
  }
  int resultRows = transposed ? columnCount() : rowCount();

  npy_intp dims[2] = {resultRows,ncols};
  PyObject* y_inout = PyArray_ZEROS(2,dims,typenum,true);

  // Create Maps to the underlying data
  
  Eigen::Map<Matrix<ValueType>> x_mat(reinterpret_cast<ValueType*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(x_f))),nrows,ncols);

  Eigen::Map<Matrix<ValueType>> y_mat(reinterpret_cast<ValueType*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(y_inout))),resultRows,ncols);

  for (int i = 0; i < ncols; ++i)
      this->apply(trans,
              Eigen::Ref<Vector<ValueType>>(x_mat.col(i)),
              Eigen::Ref<Vector<ValueType>>(y_mat.col(i)),
              1.0,0.0);

  Py_DECREF(x_f);
  return y_inout;
}


template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
DiscreteBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
    double eps, int maximumRank, bool interleave) const {
#ifdef WITH_AHMED
  throw std::runtime_error("DiscreteBoundaryOperator::"
                           "asDiscreteAcaBoundaryOperator(): "
                           "not implemented for operators of class " +
                           std::string(typeid(*this).name()) + ".");
#else
  throw std::runtime_error("DiscreteBoundaryOperator::"
                           "asDiscreteAcaBoundaryOperator(): "
                           "ACA operators are not supported because BEM++ "
                           "has been compiled without AHMED.");
#endif
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::dump() const {
  std::cout << asMatrix() << std::endl;
}



template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
operator+(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return op;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
operator-(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return static_cast<ValueType>(-1.) * op;
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
operator+(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op1,
          const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op2) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new DiscreteBoundaryOperatorSum<ValueType>(op1, op2));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
sum(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op1,
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op2) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new DiscreteBoundaryOperatorSum<ValueType>(op1, op2));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
operator-(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op1,
          const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op2) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new DiscreteBoundaryOperatorSum<ValueType>(
          op1, static_cast<ValueType>(-1.) * op2));
}

template <typename ValueType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>,
                        std::complex<double>>,
        ScalarType>,
    shared_ptr<DiscreteBoundaryOperator<ValueType>>>::type
operator*(ScalarType scalar,
          const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new ScaledDiscreteBoundaryOperator<ValueType>(
          static_cast<ValueType>(scalar), op));
}

template <typename ValueType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>,
                        std::complex<double>>,
        ScalarType>,
    shared_ptr<DiscreteBoundaryOperator<ValueType>>>::type
operator*(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op,
          ScalarType scalar) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new ScaledDiscreteBoundaryOperator<ValueType>(
          static_cast<ValueType>(scalar), op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
operator*(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op1,
          const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op2) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new DiscreteBoundaryOperatorComposition<ValueType>(op1, op2));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
mul(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op1,
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op2) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new DiscreteBoundaryOperatorComposition<ValueType>(op1, op2));
}

template <typename ValueType, typename ScalarType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
operator/(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op,
          ScalarType scalar) {
  if (scalar == static_cast<ScalarType>(0.))
    throw std::runtime_error("operator/(DiscreteBoundaryOperator, scalar): "
                             "Division by zero");
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new ScaledDiscreteBoundaryOperator<ValueType>(
          static_cast<ValueType>(static_cast<ScalarType>(1.) / scalar), op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
transpose(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new TransposedDiscreteBoundaryOperator<ValueType>(TRANSPOSE, op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
conjugate(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new TransposedDiscreteBoundaryOperator<ValueType>(CONJUGATE, op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>> conjugateTranspose(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new TransposedDiscreteBoundaryOperator<ValueType>(CONJUGATE_TRANSPOSE,
                                                        op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType>>
transpose(TranspositionMode trans,
          const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  return shared_ptr<DiscreteBoundaryOperator<ValueType>>(
      new TransposedDiscreteBoundaryOperator<ValueType>(trans, op));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(VALUE)                                      \
  template shared_ptr<const DiscreteBoundaryOperator<VALUE>> operator+(        \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);            \
  template shared_ptr<const DiscreteBoundaryOperator<VALUE>> operator-(        \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);            \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> operator+(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op1,            \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op2);           \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> sum(                    \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op1,            \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op2);           \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> operator-(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op1,            \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op2);           \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> operator*(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op1,            \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op2);           \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> mul(                    \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op1,            \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op2);           \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> transpose(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);            \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> conjugate(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);            \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> conjugateTranspose(     \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);            \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> transpose(              \
      TranspositionMode trans,                                                 \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);

#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(VALUE, SCALAR)                  \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> operator*(              \
      SCALAR scalar,                                                           \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op);            \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> operator*(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op,             \
      SCALAR scalar);                                                          \
  template shared_ptr<DiscreteBoundaryOperator<VALUE>> operator/(              \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &op,             \
      SCALAR scalar);

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, double);
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>,
                                       std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>,
                                       std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_KERNELS) ||                                        \
     defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>,
                                       std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>,
                                       std::complex<double>);
#endif

} // namespace Bempp
