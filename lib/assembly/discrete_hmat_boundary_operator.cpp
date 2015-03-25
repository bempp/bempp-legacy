// Copyright (C) 2011-2014 by the BEM++ Authors
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

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "discrete_hmat_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include <boost/numeric/conversion/converter.hpp>
#include "../hmat/compressed_matrix.hpp"
#include "../hmat/hmatrix.hpp"



namespace Bempp {

template <typename ValueType>
DiscreteHMatBoundaryOperator<ValueType>::DiscreteHMatBoundaryOperator(
    const shared_ptr<hmat::DefaultHMatrixType<ValueType>> &hMatrix)
    : m_hMatrix(hMatrix),
      m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(
          hMatrix->columns())),
      m_rangeSpace(
          Thyra::defaultSpmdVectorSpace<ValueType>(hMatrix->rows())) {}

template <typename ValueType>
unsigned int DiscreteHMatBoundaryOperator<ValueType>::rowCount() const {

  return boost::numeric::converter<unsigned int, std::size_t>::convert(
      m_hMatrix->rows());
}

template <typename ValueType>
unsigned int DiscreteHMatBoundaryOperator<ValueType>::columnCount() const {

  return boost::numeric::converter<unsigned int, std::size_t>::convert(
      m_hMatrix->columns());
}

template <typename ValueType>
shared_ptr<const hmat::DefaultHMatrixType<ValueType>> 
DiscreteHMatBoundaryOperator<ValueType>::hMatrix() const
{
    return m_hMatrix;
}

template <typename ValueType>
void DiscreteHMatBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {}

template <typename ValueType>
void DiscreteHMatBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {

  hmat::TransposeMode hmatTrans;
  if (trans == TranspositionMode::NO_TRANSPOSE)
    hmatTrans = hmat::NOTRANS;
  else if (trans == TranspositionMode::TRANSPOSE)
    hmatTrans = hmat::TRANS;
  else if (trans == TranspositionMode::CONJUGATE)
    hmatTrans = hmat::CONJ;
  else
    hmatTrans = hmat::CONJTRANS;
  Matrix<ValueType> x_inMat = x_in;
  Matrix<ValueType> y_inoutMat = y_inout;

  m_hMatrix->apply(x_inMat, y_inoutMat, hmatTrans, alpha, beta);

  y_inout = y_inoutMat.col(0);
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteHMatBoundaryOperator<ValueType>::domain() const {
  return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteHMatBoundaryOperator<ValueType>::range() const {
  return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteHMatBoundaryOperator<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS ||
          M_trans == Thyra::CONJTRANS);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteHMatBoundaryOperator);
}


