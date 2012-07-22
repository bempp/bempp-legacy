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

#include "discrete_blocked_boundary_operator.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <numeric>
#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif // WITH_TRILINOS

namespace Bempp
{

template <typename ValueType>
DiscreteBlockedBoundaryOperator<ValueType>::
DiscreteBlockedBoundaryOperator(
        const Fiber::_2dArray<shared_ptr<const Base> >& blocks,
        const std::vector<size_t>& rowCounts,
        const std::vector<size_t>& columnCounts) :
    m_blocks(blocks), m_rowCounts(rowCounts), m_columnCounts(columnCounts)
{    
    if (rowCounts.size() != blocks.extent(0))
        throw std::invalid_argument(
                "DiscreteBlockedBoundaryOperator::DiscreteBlockedBoundaryOperator(): "
                "rowCounts has incorrect length");
    if (columnCounts.size() != blocks.extent(1))
        throw std::invalid_argument(
                "DiscreteBlockedBoundaryOperator::DiscreteBlockedBoundaryOperator(): "
                "columnCounts has incorrect length");

    for (size_t col = 0; col < blocks.extent(1); ++col)
        for (size_t row = 0; row < blocks.extent(0); ++row)
            if (blocks(row, col)) {
                if (blocks(row, col)->rowCount() != rowCounts[row] ||
                        blocks(row, col)->columnCount() != columnCounts[col])
                    throw std::invalid_argument(
                            "DiscreteBlockedBoundaryOperator::"
                            "DiscreteBlockedBoundaryOperator(): "
                            "block (" + toString(row) + ", " + toString(col) + ") "
                            "has incorrect dimensions");
            }
#ifdef WITH_TRILINOS
    m_domainSpace = Thyra::defaultSpmdVectorSpace<ValueType>(
                std::accumulate(m_columnCounts.begin(), m_columnCounts.end(), 0));
    m_rangeSpace = Thyra::defaultSpmdVectorSpace<ValueType>(
                std::accumulate(m_rowCounts.begin(), m_rowCounts.end(), 0));
#endif
}

template <typename ValueType>
unsigned int
DiscreteBlockedBoundaryOperator<ValueType>::rowCount() const
{
    return std::accumulate(m_rowCounts.begin(), m_rowCounts.end(), 0);
}

template <typename ValueType>
unsigned int
DiscreteBlockedBoundaryOperator<ValueType>::columnCount() const
{
    return std::accumulate(m_columnCounts.begin(), m_columnCounts.end(), 0);
}

template <typename ValueType>
void DiscreteBlockedBoundaryOperator<ValueType>::addBlock(
        const std::vector<int>& rows,
        const std::vector<int>& cols,
        const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    throw std::runtime_error(
                "DiscreteBlockedBoundaryOperator::DiscreteBlockedBoundaryOperator(): "
                "addBlock: not implemented yet");
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteBlockedBoundaryOperator<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteBlockedBoundaryOperator<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteBlockedBoundaryOperator<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    for (size_t col = 0; col < m_blocks.extent(1); ++col)
        for (size_t row = 0; row < m_blocks.extent(0); ++row)
            if (m_blocks(row, col) && !m_blocks(row, col)->opSupported(M_trans))
                return false;
    return true;
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteBlockedBoundaryOperator<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
    bool transpose = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);
    for (int row = 0, y_start = 0; row < m_rowCounts.size(); ++row) {
        size_t y_chunk_size = m_rowCounts[row];
        arma::Col<ValueType> y_chunk(&y_inout[y_start], y_chunk_size,
                                     false /* copy_aux_mem */);
        for (int col = 0, x_start = 0; col < m_columnCounts.size(); ++col) {
            size_t x_chunk_size = m_columnCounts[col];
            shared_ptr<const Base> op =
                    transpose ? m_blocks(col, row) : m_blocks(row, col);
//            arma::Col<ValueType> x_chunk(&x_in[colStart], m_columnCounts[col],
//                                         false /* copy_aux_mem */);
            if (col == 0) {
                // This branch ensures that the "y += beta * y" part is done
                if (op)
//                    op->apply(trans, x_chunk, y_chunk, alpha, beta);
                    op->apply(trans, x_in.rows(x_start, x_start + x_chunk_size - 1),
                              y_chunk, alpha, beta);
                else {
                    if (beta == static_cast<ValueType>(0.))
                        y_chunk.fill(0.);
                    else
                        y_chunk *= beta;
                }
            }
            else
                if (op)
//                    op->apply(trans, x_chunk, y_chunk, alpha, 1.);
                    op->apply(trans, x_in.rows(x_start, x_start + x_chunk_size - 1),
                              y_chunk, alpha, 1.);
            x_start += x_chunk_size;
        }
        y_start += y_chunk_size;
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBlockedBoundaryOperator);

} // namespace Bempp
