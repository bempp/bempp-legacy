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

#include <iostream>

#include "discrete_boundary_operator.hpp"

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
void DiscreteBoundaryOperator<ValueType>::apply(const TranspositionMode trans,
                                                const Matrix<ValueType> &x_in,
                                                Matrix<ValueType> &y_inout,
                                                const ValueType alpha,
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

    applyBuiltInImpl(trans, Eigen::Ref<Vector<ValueType>>(
                                const_cast<Matrix<ValueType> &>(x_in).col(i)),
                     Eigen::Ref<Vector<ValueType>>(
                         const_cast<Matrix<ValueType> &>(y_inout).col(i)),
                     alpha, beta);
  }
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::apply(const TranspositionMode trans,
                                                const Vector<ValueType> &x_in,
                                                Vector<ValueType> &y_inout,
                                                const ValueType alpha,
                                                const ValueType beta) const {

  this->apply(
      trans,
      Eigen::Ref<Vector<ValueType>>(const_cast<Vector<ValueType> &>(x_in)),
      Eigen::Ref<Vector<ValueType>>(const_cast<Vector<ValueType> &>(y_inout)),
      alpha, beta);
}

template <typename ValueType>
Matrix<ValueType> DiscreteBoundaryOperator<ValueType>::apply(
    const TranspositionMode trans, const Matrix<ValueType> &x_in) const {

  bool transposed = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);

  unsigned int rows = (transposed ? columnCount() : rowCount());
  unsigned int cols = (transposed ? rowCount() : columnCount());

  Matrix<ValueType> result(rows, x_in.cols());

  apply(trans, x_in, result, 1.0, 0.0);
  return result;
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::apply(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {

  applyBuiltInImpl(trans, x_in, y_inout, alpha, beta);
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::dump() const {
  std::cout << asMatrix() << std::endl;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBoundaryOperator);

} // namespace Bempp
