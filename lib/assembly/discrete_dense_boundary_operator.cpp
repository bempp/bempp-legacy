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

#include "discrete_dense_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "fiber/scalar_traits.hpp"

#include <iostream>
#include <stdexcept>
#include <array>
#include <numpy/npy_common.h>

namespace Bempp {

template <typename ValueType>
DiscreteDenseBoundaryOperator<ValueType>::DiscreteDenseBoundaryOperator(
    const Matrix<ValueType> &mat)
    : m_mat(mat) {}

template <typename ValueType>
DiscreteDenseBoundaryOperator<ValueType>::DiscreteDenseBoundaryOperator() {}

template <typename ValueType>
void DiscreteDenseBoundaryOperator<ValueType>::dump() const {
  std::cout << m_mat << std::endl;
}

template <typename ValueType>
Matrix<ValueType> DiscreteDenseBoundaryOperator<ValueType>::asMatrix() const {
  return m_mat;
}

template <typename ValueType>
Matrix<ValueType>& DiscreteDenseBoundaryOperator<ValueType>::matrix() const {
  return m_mat;
}

template <typename ValueType>
unsigned int DiscreteDenseBoundaryOperator<ValueType>::rowCount() const {
  return m_mat.rows();
}

template <typename ValueType>
unsigned int DiscreteDenseBoundaryOperator<ValueType>::columnCount() const {
  return m_mat.cols();
}

template <typename ValueType>
void DiscreteDenseBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  if (block.rows() != rows.size() || block.cols() != cols.size())
    throw std::invalid_argument("DiscreteDenseBoundaryOperator::addBlock(): "
                                "incorrect block size");
  for (size_t col = 0; col < cols.size(); ++col)
    for (size_t row = 0; row < rows.size(); ++row)
      block(row, col) += alpha * m_mat(rows[row], cols[col]);
}

template <typename ValueType>
PyObject *DiscreteDenseBoundaryOperator<ValueType>::asNumpyObject() const {

  int nd = 2;
  std::array<npy_intp, 2> dims{{this->rowCount(), this->columnCount()}};

  PyGILState_STATE gstate;
  gstate = PyGILState_Ensure();
  PyObject *out =
      PyArray_New(&PyArray_Type, nd, dims.data(),
                  Fiber::ScalarTraits<ValueType>::NumpyTypeNum, NULL,
                  m_mat.data(), 0, NPY_ARRAY_F_CONTIGUOUS, NULL);
  PyGILState_Release(gstate);
  return out;
}

template <typename ValueType>
void DiscreteDenseBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {
  if (beta == static_cast<ValueType>(0.))
    y_inout.fill(static_cast<ValueType>(0.));
  else
    y_inout *= beta;

  switch (trans) {
  case NO_TRANSPOSE:
    y_inout += alpha * m_mat * x_in;
    break;
  case CONJUGATE:
    y_inout += alpha * (m_mat.conjugate()) * x_in;
    break;
  case TRANSPOSE:
    y_inout += alpha * m_mat.transpose() * x_in;
    break;
  case CONJUGATE_TRANSPOSE:
    y_inout += alpha * m_mat.adjoint() * x_in;
    break;
  default:
    throw std::invalid_argument(
        "DiscreteDenseBoundaryOperator::applyBuiltInImpl(): "
        "invalid transposition mode");
  }
}

template <typename ValueType>
shared_ptr<DiscreteDenseBoundaryOperator<ValueType>>
discreteDenseBoundaryOperator(const Matrix<ValueType> &mat) {
  typedef DiscreteDenseBoundaryOperator<ValueType> Op;
  return boost::make_shared<Op>(mat);
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(VALUE)                               \
  template shared_ptr<DiscreteDenseBoundaryOperator<VALUE>>                    \
  discreteDenseBoundaryOperator(const Matrix<VALUE> &)
FIBER_ITERATE_OVER_VALUE_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteDenseBoundaryOperator);

} // namespace Bempp
