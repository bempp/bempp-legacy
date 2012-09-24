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

#include "discrete_dense_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <stdexcept>

#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

template <typename ValueType>
DiscreteDenseBoundaryOperator<ValueType>::
DiscreteDenseBoundaryOperator(const arma::Mat<ValueType>& mat) :
    m_mat(mat)
#ifdef WITH_TRILINOS
  , m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(mat.n_cols)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(mat.n_rows))
#endif
{
}

template <typename ValueType>
void DiscreteDenseBoundaryOperator<ValueType>::dump() const
{
    std::cout << m_mat << std::endl;
}

template <typename ValueType>
arma::Mat<ValueType> DiscreteDenseBoundaryOperator<ValueType>::asMatrix() const
{
    return m_mat;
}

template <typename ValueType>
unsigned int DiscreteDenseBoundaryOperator<ValueType>::rowCount() const
{
    return m_mat.n_rows;
}

template <typename ValueType>
unsigned int DiscreteDenseBoundaryOperator<ValueType>::columnCount() const
{
    return m_mat.n_cols;
}

template <typename ValueType>
void DiscreteDenseBoundaryOperator<ValueType>::addBlock(
        const std::vector<int>& rows,
        const std::vector<int>& cols,
        const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    if (block.n_rows != rows.size() || block.n_cols != cols.size())
        throw std::invalid_argument(
                "DiscreteDenseBoundaryOperator::addBlock(): "
                "incorrect block size");
    for (size_t col = 0; col < cols.size(); ++col)
        for (size_t row = 0; row < rows.size(); ++row)
            block(row, col) += alpha*m_mat(rows[row], cols[col]);
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteDenseBoundaryOperator<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteDenseBoundaryOperator<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteDenseBoundaryOperator<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS
            || M_trans == Thyra::CONJ || M_trans == Thyra::CONJTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteDenseBoundaryOperator<ValueType>::applyBuiltInImpl(
        const TranspositionMode trans,
        const arma::Col<ValueType>& x_in,
        arma::Col<ValueType>& y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    if (beta == static_cast<ValueType>(0.))
        y_inout.fill(static_cast<ValueType>(0.));
    else
        y_inout *= beta;

    switch (trans)
    {
    case NO_TRANSPOSE:
        y_inout += alpha * m_mat * x_in;
        break;
    case CONJUGATE:
        y_inout += alpha * arma::conj(m_mat) * x_in;
        break;
    case TRANSPOSE:
        y_inout += alpha * m_mat.st() * x_in;
        break;
    case CONJUGATE_TRANSPOSE:
        y_inout += alpha * m_mat.t() * x_in;
        break;
    default:
        throw std::invalid_argument(
                "AcaApproximateLuInverse::applyBuiltInImpl(): "
                "invalid transposition mode");
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteDenseBoundaryOperator);

} // namespace Bempp
