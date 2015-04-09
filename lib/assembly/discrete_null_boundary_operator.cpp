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

#include "discrete_null_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../common/eigen_support.hpp"

#include <iostream>
#include <stdexcept>

#ifdef WITH_TRILINOS
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>
#endif

namespace Bempp {

template <typename ValueType>
DiscreteNullBoundaryOperator<ValueType>::DiscreteNullBoundaryOperator(
    size_t rowCount_, size_t columnCount_)
    :
#ifdef WITH_TRILINOS
      m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount_)),
      m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount_))
#else
      m_rowCount(rowCount_),
      m_columnCount(columnCount_)
#endif
{
}

template <typename ValueType>
void DiscreteNullBoundaryOperator<ValueType>::dump() const {
  std::cout << asMatrix() << std::endl;
}

template <typename ValueType>
Matrix<ValueType> DiscreteNullBoundaryOperator<ValueType>::asMatrix() const {
  return Matrix<ValueType>::Zero(rowCount(), columnCount());
}

template <typename ValueType>
unsigned int DiscreteNullBoundaryOperator<ValueType>::rowCount() const {
#ifdef WITH_TRILINOS
  return m_rangeSpace->dim();
#else
  return m_rowCount;
#endif
}

template <typename ValueType>
unsigned int DiscreteNullBoundaryOperator<ValueType>::columnCount() const {
#ifdef WITH_TRILINOS
  return m_domainSpace->dim();
#else
  return m_columnCount;
#endif
}

template <typename ValueType>
void DiscreteNullBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  // don't do anything
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteNullBoundaryOperator<ValueType>::domain() const {
  return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteNullBoundaryOperator<ValueType>::range() const {
  return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteNullBoundaryOperator<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS ||
          M_trans == Thyra::CONJ || M_trans == Thyra::CONJTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteNullBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {
  if (beta == static_cast<ValueType>(0.))
    y_inout.fill(static_cast<ValueType>(0.));
  else
    y_inout *= beta;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteNullBoundaryOperator);

} // namespace Bempp
