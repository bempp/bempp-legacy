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

#include "bempp/common/config_trilinos.hpp"

#include "discrete_boundary_operator_sum.hpp"
#include "discrete_aca_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename ValueType>
DiscreteBoundaryOperatorSum<ValueType>::DiscreteBoundaryOperatorSum(
    const shared_ptr<const Base> &term1, const shared_ptr<const Base> &term2)
    : m_term1(term1), m_term2(term2) {
  if (!m_term1 || !m_term2)
    throw std::invalid_argument(
        "DiscreteBoundaryOperatorSum::DiscreteBoundaryOperatorSum(): "
        "arguments must not be NULL");
  if (m_term1->rowCount() != m_term2->rowCount() ||
      m_term1->columnCount() != m_term2->columnCount())
    throw std::invalid_argument(
        "DiscreteBoundaryOperatorSum::DiscreteBoundaryOperatorSum(): "
        "both terms must have the same dimensions");
  // TODO: perhaps test for compatibility of Thyra spaces
}

template <typename ValueType>
arma::Mat<ValueType> DiscreteBoundaryOperatorSum<ValueType>::asMatrix() const {
  arma::Mat<ValueType> result = m_term1->asMatrix();
  result += m_term2->asMatrix();
  return result;
}

template <typename ValueType>
unsigned int DiscreteBoundaryOperatorSum<ValueType>::rowCount() const {
  return m_term1->rowCount();
}

template <typename ValueType>
unsigned int DiscreteBoundaryOperatorSum<ValueType>::columnCount() const {
  return m_term1->columnCount();
}

template <typename ValueType>
void DiscreteBoundaryOperatorSum<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, arma::Mat<ValueType> &block) const {
  m_term1->addBlock(rows, cols, alpha, block);
  m_term2->addBlock(rows, cols, alpha, block);
}

#ifdef WITH_AHMED
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
DiscreteBoundaryOperatorSum<ValueType>::asDiscreteAcaBoundaryOperator(
    double eps, int maximumRank, bool interleave) const {
  typedef DiscreteAcaBoundaryOperator<ValueType> DiscreteAcaOp;
  shared_ptr<const DiscreteAcaOp> acaOp1 =
      boost::dynamic_pointer_cast<const DiscreteAcaOp>(
          m_term1->asDiscreteAcaBoundaryOperator(eps, maximumRank, interleave));
  shared_ptr<const DiscreteAcaOp> acaOp2 =
      boost::dynamic_pointer_cast<const DiscreteAcaOp>(
          m_term2->asDiscreteAcaBoundaryOperator(eps, maximumRank, interleave));

  if (maximumRank < 0)
    maximumRank = std::max(acaOp1->maximumRank(), acaOp2->maximumRank());
  if (eps < 0)
    eps = std::min(acaOp1->eps(), acaOp2->eps());

  return acaOperatorSum<ValueType>(acaOp1, acaOp2, eps, maximumRank);
}
#endif // WITH_AHMED

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteBoundaryOperatorSum<ValueType>::domain() const {
  return m_term1->domain();
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteBoundaryOperatorSum<ValueType>::range() const {
  return m_term1->range();
}

template <typename ValueType>
bool DiscreteBoundaryOperatorSum<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  return (m_term1->opSupported(M_trans) && m_term2->opSupported(M_trans));
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteBoundaryOperatorSum<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const arma::Col<ValueType> &x_in,
    arma::Col<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  m_term1->apply(trans, x_in, y_inout, alpha, beta);
  m_term2->apply(trans, x_in, y_inout, alpha,
                 1. /* "+ beta * y_inout" has already been done */);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBoundaryOperatorSum);

} // namespace Bempp
