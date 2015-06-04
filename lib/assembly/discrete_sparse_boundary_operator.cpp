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

#include "discrete_sparse_boundary_operator.hpp"

#include "../common/boost_shared_array_fwd.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/parallelization_options.hpp"

#include <iostream>
#include <stdexcept>

namespace Bempp {

// Helper functions for the applyBuiltIn member function
namespace {

template <typename ValueType>
void reallyApplyBuiltInImpl(const RealSparseMatrix &mat,
                            const TranspositionMode trans,
                            const Eigen::Ref<Vector<ValueType>> &x_in,
                            Eigen::Ref<Vector<ValueType>> y_inout,
                            const ValueType alpha, const ValueType beta);

template <>
void reallyApplyBuiltInImpl<double>(const RealSparseMatrix &mat,
                                    const TranspositionMode trans,
                                    const Eigen::Ref<Vector<double>> &x_in,
                                    Eigen::Ref<Vector<double>> y_inout,
                                    const double alpha, const double beta) {
  if (beta == 0.)
    y_inout.setZero();

  if (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE) {
    assert(mat.rows() == static_cast<int>(x_in.rows()));
    assert(mat.cols() == static_cast<int>(y_inout.rows()));
  } else {
    assert(mat.cols() == static_cast<int>(x_in.rows()));
    assert(mat.rows() == static_cast<int>(y_inout.rows()));
  }

  if (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE)
    y_inout = alpha * (mat.transpose() * x_in) + beta * y_inout;
  else
    y_inout = alpha * (mat * x_in) + beta * y_inout;
}

template <>
void reallyApplyBuiltInImpl<float>(const RealSparseMatrix &mat,
                                   const TranspositionMode trans,
                                   const Eigen::Ref<Vector<float>> &x_in,
                                   Eigen::Ref<Vector<float>> y_inout,
                                   const float alpha, const float beta) {
  // Copy the float vectors to double vectors
  Vector<double> x_in_double = x_in.cast<double>();

  Vector<double> y_inout_double = y_inout.cast<double>();

  // Do the operation on the double vectors
  reallyApplyBuiltInImpl<double>(
      mat, trans, Eigen::Ref<Vector<double>>(x_in_double),
      Eigen::Ref<Vector<double>>(y_inout_double), alpha, beta);

  // Copy the result back to the float vector

  y_inout = y_inout_double.cast<float>();
}

template <>
void reallyApplyBuiltInImpl<std::complex<float>>(
    const RealSparseMatrix &mat, const TranspositionMode trans,
    const Eigen::Ref<Vector<std::complex<float>>> &x_in,
    Eigen::Ref<Vector<std::complex<float>>> y_inout,
    const std::complex<float> alpha, const std::complex<float> beta) {
  // Do the y_inout *= beta part
  const std::complex<float> zero(0.f, 0.f);
  if (beta == zero)
    y_inout.setZero();
  else
    y_inout *= beta;

  // Separate the real and imaginary components and store them in
  // double-precision vectors

  Vector<double> x_real = x_in.real().cast<double>();
  Vector<double> x_imag = x_in.imag().cast<double>();
  Vector<double> y_real = y_inout.real().cast<double>();
  Vector<double> y_imag = y_inout.imag().cast<double>();

  //  Vector<double> x_real(x_in.rows());
  //  for (size_t i = 0; i < x_in.rows(); ++i)
  //    x_real(i) = x_in(i).real();
  //  Vector<double> x_imag(x_in.rows());
  //  for (size_t i = 0; i < x_in.rows(); ++i)
  //    x_imag(i) = x_in(i).imag();
  //  Vector<double> y_real(y_inout.rows());
  //  for (size_t i = 0; i < y_inout.rows(); ++i)
  //    y_real(i) = y_inout(i).real();
  //  Vector<double> y_imag(y_inout.rows());
  //  for (size_t i = 0; i < y_inout.rows(); ++i)
  //    y_imag(i) = y_inout(i).imag();

  // Do the "+= alpha A x" part (in steps)
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real),
                                 Eigen::Ref<Vector<double>>(y_real),
                                 alpha.real(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag),
                                 Eigen::Ref<Vector<double>>(y_real),
                                 -alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real),
                                 Eigen::Ref<Vector<double>>(y_imag),
                                 alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag),
                                 Eigen::Ref<Vector<double>>(y_imag),
                                 alpha.real(), 1.);

  y_inout.real() = y_real.cast<float>();
  y_inout.imag() = y_imag.cast<float>();

  //  // Copy the result back to the complex vector
  //  for (size_t i = 0; i < y_inout.rows(); ++i)
  //    y_inout(i) = std::complex<float>(y_real(i), y_imag(i));
}

template <>
void reallyApplyBuiltInImpl<std::complex<double>>(
    const RealSparseMatrix &mat, const TranspositionMode trans,
    const Eigen::Ref<Vector<std::complex<double>>> &x_in,
    Eigen::Ref<Vector<std::complex<double>>> y_inout,
    const std::complex<double> alpha, const std::complex<double> beta) {
  // Do the y_inout *= beta part
  const std::complex<double> zero(0.f, 0.f);
  if (beta == zero)
    y_inout.setZero();
  else
    y_inout *= beta;

  // Separate the real and imaginary components and store them in
  // double-precision vectors

  Vector<double> x_real = x_in.real();
  Vector<double> x_imag = x_in.imag();
  Vector<double> y_real = y_inout.real();
  Vector<double> y_imag = y_inout.imag();

  //  Vector<double> x_real(x_in.rows());
  //  for (size_t i = 0; i < x_in.rows(); ++i)
  //    x_real(i) = x_in(i).real();
  //  Vector<double> x_imag(x_in.rows());
  //  for (size_t i = 0; i < x_in.rows(); ++i)
  //    x_imag(i) = x_in(i).imag();
  //  Vector<double> y_real(y_inout.rows());
  //  for (size_t i = 0; i < y_inout.rows(); ++i)
  //    y_real(i) = y_inout(i).real();
  //  Vector<double> y_imag(y_inout.rows());
  //  for (size_t i = 0; i < y_inout.rows(); ++i)
  //    y_imag(i) = y_inout(i).imag();

  // Do the "+= alpha A x" part (in steps)
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real),
                                 Eigen::Ref<Vector<double>>(y_real),
                                 alpha.real(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag),
                                 Eigen::Ref<Vector<double>>(y_real),
                                 -alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real),
                                 Eigen::Ref<Vector<double>>(y_imag),
                                 alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag),
                                 Eigen::Ref<Vector<double>>(y_imag),
                                 alpha.real(), 1.);

  y_inout.real() = y_real;
  y_inout.imag() = y_imag;

  //  // Copy the result back to the complex vector
  //  for (size_t i = 0; i < y_inout.rows(); ++i)
  //    y_inout(i) = std::complex<float>(y_real(i), y_imag(i));
}

} // namespace

template <typename ValueType>
DiscreteSparseBoundaryOperator<ValueType>::DiscreteSparseBoundaryOperator(
    const shared_ptr<const RealSparseMatrix> &mat, int symmetry,
    TranspositionMode trans)
    : m_mat(mat), m_symmetry(symmetry), m_trans(trans) {}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::dump() const {
  if (isTransposed())
    std::cout << "Transpose of " << *m_mat << std::endl;
  else
    std::cout << *m_mat << std::endl;
}

template <typename ValueType>
Matrix<ValueType> DiscreteSparseBoundaryOperator<ValueType>::asMatrix() const {

  bool transposed = isTransposed();
  if (transposed)
    return m_mat->transpose().template cast<ValueType>();
  else
    return m_mat->cast<ValueType>();
}

template <typename ValueType>
unsigned int DiscreteSparseBoundaryOperator<ValueType>::rowCount() const {
  return isTransposed() ? m_mat->cols() : m_mat->rows();
}

template <typename ValueType>
unsigned int DiscreteSparseBoundaryOperator<ValueType>::columnCount() const {
  return isTransposed() ? m_mat->rows() : m_mat->cols();
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {

  throw std::runtime_error(
      "DiscreteSparseBoundaryOperator::addBlock(): not implemented.");
}

template <typename ValueType>
shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>>
DiscreteSparseBoundaryOperator<ValueType>::castToSparse(const shared_ptr<
    const DiscreteBoundaryOperator<ValueType>> &discreteOperator) {
  shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>> result =
      boost::dynamic_pointer_cast<
          const DiscreteSparseBoundaryOperator<ValueType>>(discreteOperator);
  if (result.get() == 0 && (discreteOperator.get() != 0))
    throw std::bad_cast();
  return result;
}

template <typename ValueType>
shared_ptr<const RealSparseMatrix>
DiscreteSparseBoundaryOperator<ValueType>::sparseMatrix() const {
  return m_mat;
}

template <typename ValueType>
TranspositionMode
DiscreteSparseBoundaryOperator<ValueType>::transpositionMode() const {
  return m_trans;
}

template <typename ValueType>
bool DiscreteSparseBoundaryOperator<ValueType>::isTransposed() const {
  return m_trans & (TRANSPOSE | CONJUGATE_TRANSPOSE);
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {
  TranspositionMode realTrans = trans;
  bool transposed = isTransposed();
  if (transposed)
    switch (trans) {
    case NO_TRANSPOSE:
      realTrans = TRANSPOSE;
      break;
    case TRANSPOSE:
      realTrans = NO_TRANSPOSE;
      break;
    case CONJUGATE:
      realTrans = CONJUGATE_TRANSPOSE;
      break;
    case CONJUGATE_TRANSPOSE:
      realTrans = CONJUGATE;
      break;
      // default: should not happen; anyway, don't change trans
    }

  reallyApplyBuiltInImpl(*m_mat, realTrans, x_in, y_inout, alpha, beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteSparseBoundaryOperator);

} // namespace Bempp
