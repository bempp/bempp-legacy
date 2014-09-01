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

#include "discrete_hmat_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include <boost/numeric/conversion/converter.hpp>

namespace Bempp {

template <typename ValueType>
DiscreteHMatBoundaryOperator<ValueType>::DiscreteHMatBoundaryOperator(
    const shared_ptr<hmat::CompressedMatrix<ValueType>> &compressedMatrix)
    : m_compressedMatrix(compressedMatrix),
      m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(
          compressedMatrix->columns())),
      m_rangeSpace(
          Thyra::defaultSpmdVectorSpace<ValueType>(compressedMatrix->rows())) {}

template <typename ValueType>
unsigned int DiscreteHMatBoundaryOperator<ValueType>::rowCount() const {

  return boost::numeric::converter<unsigned int, std::size_t>::convert(
      m_compressedMatrix->rows());
}

template <typename ValueType>
unsigned int DiscreteHMatBoundaryOperator<ValueType>::columnCount() const {

  return boost::numeric::converter<unsigned int, std::size_t>::convert(
      m_compressedMatrix->columns());
}

template <typename ValueType>
void DiscreteHMatBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, arma::Mat<ValueType> &block) const {}

template <typename ValueType>
void DiscreteHMatBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const arma::Col<ValueType> &x_in,
    arma::Col<ValueType> &y_inout, const ValueType alpha,
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
  m_compressedMatrix->apply(x_in, y_inout, hmatTrans, alpha, beta);
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
