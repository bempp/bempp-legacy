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
#include "bempp/common/config_ahmed.hpp"

#include "discrete_blocked_boundary_operator.hpp"

#include "discrete_aca_boundary_operator.hpp"
#include "ahmed_aux.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <numeric>
#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif // WITH_TRILINOS

namespace Bempp
{

namespace
{

template <typename ValueType>
void copySonsAdjustingIndices(
        const typename DiscreteAcaBoundaryOperator<ValueType>::
            AhmedBemBlcluster* source,
        typename DiscreteAcaBoundaryOperator<ValueType>::
            AhmedBemBlcluster* dest,
        size_t rowOffset,
        size_t columnOffset,
        size_t indexOffset)
{
    typedef typename DiscreteAcaBoundaryOperator<ValueType>::
            AhmedBemBlcluster AhmedBemBlcluster;

    unsigned int sonCount = source->getns();
    if (sonCount > 0) {
        // how I miss auto_array...
        std::vector<blcluster*> destSons(sonCount, 0);
        try {
            for (unsigned int row = 0, i = 0; row < source->getnrs(); ++row)
                for (unsigned int col = 0; col < source->getncs(); ++col, ++i) {
                    if (const AhmedBemBlcluster* sourceSon =
                            dynamic_cast<const AhmedBemBlcluster*>(
                                source->getson(row, col))) {
                        std::auto_ptr<AhmedBemBlcluster> destSon(
                                    new AhmedBemBlcluster(
                                        rowOffset + sourceSon->getb1(),
                                        columnOffset + sourceSon->getb2(),
                                        sourceSon->getn1(),
                                        sourceSon->getn2()));
                        copySonsAdjustingIndices<ValueType>(
                                    sourceSon, destSon.get(),
                                    rowOffset, columnOffset, indexOffset);
                        destSons[i] = destSon.release();
                    }
                }
            dest->setsons(source->getnrs(), source->getncs(), &destSons[0]);
        }
        catch (...) {
            for (unsigned int i = 0; i < sonCount; ++i)
                delete destSons[i];
        }
    } else {
        dest->setidx(indexOffset + source->getidx());
        dest->setadm(dest->isadm());
        dest->setsep(dest->issep());
    }
}

} // namespace

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
shared_ptr<const DiscreteBoundaryOperator<ValueType> >
DiscreteBlockedBoundaryOperator<ValueType>::getComponent(int row, int col) const{
    return m_blocks(row,col);
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

template <typename ValueType>
shared_ptr<const DiscreteBlockedBoundaryOperator<ValueType> >
DiscreteBlockedBoundaryOperator<ValueType>::asDiscreteAcaBlockedBoundaryOperator(
        double eps, int maximumRank) const
{
#ifdef WITH_AHMED
    size_t nrows = m_blocks.extent(0);
    size_t ncols = m_blocks.extent(1);
    Fiber::_2dArray<shared_ptr<const Base> > acaBlocks(nrows, ncols);
    for (size_t i = 0; i < nrows; i++){
        for (size_t j = 0; j < ncols; j++){
            if (m_blocks(i, j).get() != 0)
                acaBlocks(i, j) =
                    m_blocks(i, j)->asDiscreteAcaBoundaryOperator(eps, maximumRank);
        }
    }
    return shared_ptr<const DiscreteBlockedBoundaryOperator<ValueType> >(
                new DiscreteBlockedBoundaryOperator(
                    acaBlocks, m_rowCounts, m_columnCounts));
#else // WITH_AHMED
    throw std::runtime_error("DiscreteBlockedBoundaryOperator::"
                             "asDiscreteAcaBlockedBoundaryOperator(): "
                             "ACA operators are not supported because BEM++ "
                             "has been compiled without AHMED.");
#endif
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> >
DiscreteBlockedBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
        double eps, int maximumRank) const
{
#ifdef WITH_AHMED
    // Get block counts
    const size_t blockRowCount = m_blocks.extent(0);
    const size_t blockColCount = m_blocks.extent(1);

    // Convert blocks into ACA operators
    typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
    typedef typename AcaOp::AhmedBemBlcluster AhmedBemBlcluster;
    typedef typename AcaOp::AhmedMblock AhmedMblock;
    typedef typename AcaOp::AhmedMblockArray AhmedMblockArray;
    typedef typename AcaOp::AhmedConstMblockArray AhmedConstMblockArray;

    Fiber::_2dArray<shared_ptr<const AcaOp> > acaBlocks(blockRowCount,
                                                        blockColCount);
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (m_blocks(row, col))
                acaBlocks(row, col) = boost::shared_dynamic_cast<const AcaOp>(
                            m_blocks(row, col)->asDiscreteAcaBoundaryOperator(
                                eps, maximumRank));
            else
                throw std::runtime_error("DiscreteBlockedBoundaryOperator::"
                                         "asDiscreteAcaBoundaryOperator(): "
                                         "operators with empty blocks are not "
                                         "supported yet");

    // Collect shared pointers to all mblock arrays (they must be stored
    // in the newly created operator so that the mblocks don't get deleted
    // prematurely). Also calculate mblock offsets.
    std::vector<AhmedConstMblockArray> sharedMblocks;
    Fiber::_2dArray<size_t> indexOffsets(blockRowCount, blockColCount);
    std::fill(indexOffsets.begin(), indexOffsets.end(), 0);
    size_t totalMblockCount = 0;
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (acaBlocks(row, col)) {
                sharedMblocks.push_back(acaBlocks(row, col)->blocks());
                sharedMblocks.insert(sharedMblocks.end(),
                                     acaBlocks(row, col)->sharedBlocks().begin(),
                                     acaBlocks(row, col)->sharedBlocks().end());
                indexOffsets(row, col) = totalMblockCount;
                totalMblockCount += acaBlocks(row, col)->blockCount();
            }

    // Create a shared array of pointers to all mblocks. Note: the array owns
    // only pointers, not mblocks themselves.
    AhmedMblockArray mblocks(new AhmedMblock*[totalMblockCount]);
    for (size_t col = 0, offset = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (acaBlocks(row, col)) {
                size_t localMblockCount = acaBlocks(row, col)->blockCount();
                std::copy(&acaBlocks(row, col)->blocks()[0],
                          &acaBlocks(row, col)->blocks()[localMblockCount],
                          &mblocks[offset]);
                offset += localMblockCount;
            }

    // Get row and column offsets and total number of rows and columns
    std::vector<size_t> rowOffsets(blockRowCount, 0);
    for (size_t i = 1; i < blockRowCount; ++i)
        rowOffsets[i] = rowOffsets[i - 1] + m_rowCounts[i - 1];
    const size_t totalRowCount = rowOffsets.back() + m_rowCounts.back();
    std::vector<size_t> colOffsets(blockColCount, 0);
    for (size_t i = 1; i < blockColCount; ++i)
        colOffsets[i] = colOffsets[i - 1] + m_columnCounts[i - 1];
    const size_t totalColCount = colOffsets.back() + m_columnCounts.back();

    // Recursively copy block cluster trees, adjusting indices
    // This vector stores a row-major array (as expected by AHMED)
    std::vector<blcluster*> roots(blockRowCount * blockColCount, 0);
    try {
        for (size_t row = 0, i = 0; row < blockRowCount; ++row)
            for (size_t col = 0; col < blockColCount; ++col, ++i)
                if (acaBlocks(row, col)) {
                    std::auto_ptr<AhmedBemBlcluster> root(
                                new AhmedBemBlcluster(
                                    rowOffsets[row], colOffsets[col],
                                    m_rowCounts[row], m_columnCounts[col]));
                    copySonsAdjustingIndices<ValueType>(
                                acaBlocks(row, col)->blockCluster(),
                                root.get(),
                                rowOffsets[row],
                                colOffsets[col],
                                indexOffsets(row, col));
                    roots[i] = root.release();
                }
    }
    catch (...) {
        for (size_t i = 0; i < roots.size(); ++i)
            delete roots[i];
        throw; // rethrow exception
    }

    // Create the root node of the combined cluster tree and set its sons
    // (TODO: see what happens if any sons are NULL)
    std::auto_ptr<AhmedBemBlcluster> mblockCluster(
                new AhmedBemBlcluster(0, 0, totalRowCount, totalColCount));
    mblockCluster->setsons(blockRowCount, blockColCount, &roots[0]);

    // Create the domain permutation array
    std::vector<unsigned int> o2pDomain;
    o2pDomain.reserve(totalRowCount);
    for (size_t col = 0; col < blockColCount; ++col) {
        int chosenRow = -1;
        for (size_t row = 0; row < blockRowCount; ++row)
            if (acaBlocks(row, col)) {
                if (chosenRow < 0)
                    chosenRow = row;
                else { // ensure that all permutations are compatible
                    if (acaBlocks(chosenRow, col)->domainPermutation() !=
                            acaBlocks(row, col)->domainPermutation())
                        throw std::runtime_error(
                                "DiscreteBlockedBoundaryOperator::"
                                "asDiscreteAcaBoundaryOperator(): "
                                "Domain index permutations of blocks (" +
                                toString(chosenRow) + ", " + toString(col) +
                                ") and (" +
                                toString(row) + ", " + toString(col) +
                                ") are not equal");
                }
            }
        std::vector<unsigned int> local_o2p =
                acaBlocks(chosenRow, col)->domainPermutation().permutedIndices();
        for (size_t i = 0; i < local_o2p.size(); ++i)
            o2pDomain.push_back(local_o2p[i] + colOffsets[col]);
    }
    assert(o2pDomain.size() == totalColCount);

    // Create the domain permutation array
    std::vector<unsigned int> o2pRange;
    o2pRange.reserve(totalColCount);
    for (size_t row = 0; row < blockRowCount; ++row) {
        int chosenCol = -1;
        for (size_t col = 0; col < blockColCount; ++col)
            if (acaBlocks(row, col)) {
                if (chosenCol < 0)
                    chosenCol = col;
                else { // ensure that all permutations are compatible
                    if (acaBlocks(row, chosenCol)->rangePermutation() !=
                            acaBlocks(row, col)->rangePermutation())
                        throw std::runtime_error(
                                "DiscreteBlockedBoundaryOperator::"
                                "asDiscreteAcaBoundaryOperator(): "
                                "Range index permutations of blocks (" +
                                toString(row) + ", " + toString(chosenCol) +
                                ") and (" +
                                toString(row) + ", " + toString(col) +
                                ") are not equal");
                }
            }
        std::vector<unsigned int> local_o2p =
                acaBlocks(row, chosenCol)->rangePermutation().permutedIndices();
        for (size_t i = 0; i < local_o2p.size(); ++i)
            o2pRange.push_back(local_o2p[i] + rowOffsets[row]);
    }
    assert(o2pRange.size() == totalRowCount);

    // Calculate the overall maximum rank
    int overallMaximumRank = 0;
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (acaBlocks(row, col))
                overallMaximumRank = std::max(overallMaximumRank,
                                       acaBlocks(row, col)->maximumRank());

    // Check symmetry
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (acaBlocks(row, col) && acaBlocks(row, col)->symmetry() != 0)
                throw std::runtime_error("DiscreteBlockedBoundaryOperator::"
                                         "asDiscreteAcaBoundaryOperator(): "
                                         "Joining of symmetric ACA operators "
                                         "is not supported yet");

    // Get some parallelization options
    Fiber::ParallelizationOptions parallelOptions;
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (acaBlocks(row, col)) {
                parallelOptions = acaBlocks(row, col)->parallelizationOptions();
            }

    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new AcaOp(
                    totalRowCount, totalColCount, overallMaximumRank,
                    0 /* symmetry */,
                    mblockCluster, mblocks,
                    o2pDomain, o2pRange,
                    parallelOptions, sharedMblocks));
    return result;

#else // WITH_AHMED
    throw std::runtime_error("DiscreteBoundaryOperator::"
                             "asDiscreteAcaBoundaryOperator(): "
                             "ACA operators are not supported because BEM++ "
                             "has been compiled without AHMED.");
#endif
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
    size_t y_count = transpose ? m_columnCounts.size() : m_rowCounts.size();
    size_t x_count = transpose ? m_rowCounts.size() : m_columnCounts.size();

    for (int yi = 0, y_start = 0; yi < y_count; ++yi) {
        size_t y_chunk_size = transpose ? m_columnCounts[yi] : m_rowCounts[yi];
        arma::Col<ValueType> y_chunk(&y_inout[y_start], y_chunk_size,
                                     false /* copy_aux_mem */);
        for (int xi = 0, x_start = 0; xi < x_count; ++xi) {
	    size_t x_chunk_size = transpose ? m_rowCounts[xi] : m_columnCounts[xi];
            shared_ptr<const Base> op =
                    transpose ? m_blocks(xi, yi) : m_blocks(yi, xi);
//            arma::Col<ValueType> x_chunk(&x_in[colStart], m_columnCounts[col],
//                                         false /* copy_aux_mem */);
            if (xi == 0) {
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
