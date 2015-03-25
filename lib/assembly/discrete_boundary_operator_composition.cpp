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

#include "discrete_boundary_operator_composition.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../common/eigen_support.hpp"

namespace Bempp {

template <typename ValueType>
DiscreteBoundaryOperatorComposition<ValueType>::
    DiscreteBoundaryOperatorComposition(const shared_ptr<const Base> &outer,
                                        const shared_ptr<const Base> &inner)
    : m_outer(outer), m_inner(inner) {
  if (!m_outer || !m_inner)
    throw std::invalid_argument("DiscreteBoundaryOperatorComposition::"
                                "DiscreteBoundaryOperatorComposition(): "
                                "arguments must not be NULL");
  if (m_outer->columnCount() != m_inner->rowCount())
    throw std::invalid_argument("DiscreteBoundaryOperatorComposition::"
                                "DiscreteBoundaryOperatorComposition(): "
                                "term dimensions do not match");
  // TODO: perhaps test for compatibility of Thyra spaces
}

template <typename ValueType>
unsigned int DiscreteBoundaryOperatorComposition<ValueType>::rowCount() const {
  return m_outer->rowCount();
}

template <typename ValueType>
unsigned int
DiscreteBoundaryOperatorComposition<ValueType>::columnCount() const {
  return m_inner->columnCount();
}

template <typename ValueType>
void DiscreteBoundaryOperatorComposition<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  throw std::runtime_error("DiscreteBoundaryOperatorComposition::"
                           "DiscreteBoundaryOperatorComposition(): "
                           "addBlock: not implemented yet");
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteBoundaryOperatorComposition<ValueType>::domain() const {
  return m_inner->domain();
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteBoundaryOperatorComposition<ValueType>::range() const {
  return m_outer->range();
}

template <typename ValueType>
bool DiscreteBoundaryOperatorComposition<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  return (m_outer->opSupported(M_trans) && m_inner->opSupported(M_trans));
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteBoundaryOperatorComposition<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
    Matrix<ValueType> x_inMat = x_in;
    Matrix<ValueType> y_inoutMat = y_inout;
  if (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE) {
    Matrix<ValueType> tmp(m_outer->columnCount(),1);
    m_outer->apply(trans, x_inMat, tmp, alpha, 0.);
    m_inner->apply(trans, tmp, y_inoutMat, 1., beta);
  } else {
    Matrix<ValueType> tmp(m_inner->rowCount(),1);
    m_inner->apply(trans, x_inMat, tmp, alpha, 0.);
    m_outer->apply(trans, tmp, y_inoutMat, 1., beta);
  }
  y_inout = y_inoutMat.col(0);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(
    DiscreteBoundaryOperatorComposition);

} // namespace Bempp
