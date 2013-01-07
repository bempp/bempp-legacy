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
#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#endif
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <numeric>
#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif // WITH_TRILINOS

namespace Bempp
{

// Functions used in
// DiscreteBlockBoundaryOperator::asDiscreteAcaBoundaryOperator().
namespace
{

#ifdef WITH_AHMED
template <typename ValueType>
void collectMblocks(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks,
        const std::vector<size_t>& rowCounts,
        const std::vector<size_t>& columnCounts,
        std::vector<typename DiscreteAcaBoundaryOperator<ValueType>::
            AhmedConstMblockArray>& allSharedMblocks,
        typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblockArray& mblocks,
        Fiber::_2dArray<size_t>& indexOffsets,
        size_t& totalMblockCount)
{
    typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
    typedef typename AcaOp::AhmedMblock AhmedMblock;
    typedef typename AcaOp::AhmedMblockArray AhmedMblockArray;
    typedef typename AcaOp::AhmedConstMblockArray AhmedConstMblockArray;

    const size_t blockRowCount = acaBlocks.extent(0);
    const size_t blockColCount = acaBlocks.extent(1);

    // Put:
    // - into explSharedMblocks: arrays of mblocks making up the H-matrix
    //   of each block (including the single-element mblock arrays created for
    //   each empty block)
    // - into allSharedMblocks: the above plus the arrays of mblocks on which
    //   the non-empty blocks depend implicitly
    // - into mblockCounts: numbers of mblocks making up the H-matrix of each
    //   block
    // - into indexOffsets: numbers by which the indices of the mblocks making
    //   up the H-matrix of each block should be offset
    // - into totalMblockCount: the total size of all arrays from
    //   explSharedBlocks

    Fiber::_2dArray<AhmedConstMblockArray> explSharedMblocks(
                blockRowCount, blockColCount);
    Fiber::_2dArray<size_t> mblockCounts(blockRowCount, blockColCount);
    indexOffsets.set_size(blockRowCount, blockColCount);
    std::fill(indexOffsets.begin(), indexOffsets.end(), 0);
    totalMblockCount = 0;
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row) {
            if (acaBlocks(row, col)) {
                const AcaOp& activeOp = *acaBlocks(row, col);
                explSharedMblocks(row, col) = activeOp.blocks();
                allSharedMblocks.insert(allSharedMblocks.end(),
                                        activeOp.sharedBlocks().begin(),
                                        activeOp.sharedBlocks().end());
                mblockCounts(row, col) = activeOp.blockCount();
            } else {
                // Allocate a single empty block
                AhmedMblockArray dummyBlocks =
                        allocateAhmedMblockArray<ValueType>(1);
                dummyBlocks[0] =
                        new AhmedMblock(rowCounts[row], columnCounts[col]);
                explSharedMblocks(row, col) = dummyBlocks;
                mblockCounts(row, col) = 1;
            }
            allSharedMblocks.push_back(explSharedMblocks(row, col));
            indexOffsets(row, col) = totalMblockCount;
            totalMblockCount += mblockCounts(row, col);
        }

    // Create a shared array of pointers to all mblocks. Note: the array owns
    // only pointers, not mblocks themselves.
    mblocks.reset(new AhmedMblock*[totalMblockCount]);
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            std::copy(&explSharedMblocks(row, col)[0],
                      &explSharedMblocks(row, col)[mblockCounts(row, col)],
                      &mblocks[indexOffsets(row, col)]);
}

template <typename ValueType>
void copySonsAdjustingIndices(
        const blcluster* source,
        blcluster* dest,
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
                    const blcluster* sourceSon = source->getson(row, col);
                    if (!sourceSon)
                        continue;
                    std::auto_ptr<blcluster> destSon;
                    if (const AhmedBemBlcluster* bemSourceSon =
                        dynamic_cast<const AhmedBemBlcluster*>(sourceSon))
                        destSon.reset(new AhmedBemBlcluster(
                                          rowOffset + sourceSon->getb1(),
                                          columnOffset + sourceSon->getb2(),
                                          sourceSon->getn1(),
                                          sourceSon->getn2()));
                    else
                        destSon.reset(new blcluster(
                                          rowOffset + sourceSon->getb1(),
                                          columnOffset + sourceSon->getb2(),
                                          sourceSon->getn1(),
                                          sourceSon->getn2()));
                    copySonsAdjustingIndices<ValueType>(
                        sourceSon, destSon.get(),
                        rowOffset, columnOffset, indexOffset);
                    destSons[i] = destSon.release();
                }
            dest->setsons(source->getnrs(), source->getncs(), &destSons[0]);
        }
        catch (...) {
            for (unsigned int i = 0; i < sonCount; ++i)
                delete destSons[i];
            throw; // rethrow
        }
    } else {
        dest->setidx(indexOffset + source->getidx());
        dest->setadm(dest->isadm());
        dest->setsep(dest->issep());
    }
}

template <typename ValueType>
double overallEps(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks)
{
    double result = std::numeric_limits<double>::max();
    for (size_t col = 0; col < acaBlocks.extent(1); ++col)
        for (size_t row = 0; row < acaBlocks.extent(0); ++row)
            if (acaBlocks(row, col))
                result = std::min(result, acaBlocks(row, col)->eps());
    return result;
}

template <typename ValueType>
int overallMaximumRank(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks)
{
    int result = 0;
    for (size_t col = 0; col < acaBlocks.extent(1); ++col)
        for (size_t row = 0; row < acaBlocks.extent(0); ++row)
            if (acaBlocks(row, col))
                result = std::max(result, acaBlocks(row, col)->maximumRank());
    return result;
}

template <typename ValueType>
int overallSymmetry(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks)
{
    for (size_t col = 0; col < acaBlocks.extent(1); ++col)
        for (size_t row = 0; row < acaBlocks.extent(0); ++row)
            if (acaBlocks(row, col) && acaBlocks(row, col)->symmetry() != 0)
                throw std::runtime_error("DiscreteBlockedBoundaryOperator::"
                                         "asDiscreteAcaBoundaryOperator(): "
                                         "Joining of symmetric ACA operators "
                                         "is not supported yet");
    return 0;
}

template <typename ValueType>
Fiber::ParallelizationOptions overallParallelizationOptions(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks)
{
    for (size_t col = 0; col < acaBlocks.extent(1); ++col)
        for (size_t row = 0; row < acaBlocks.extent(0); ++row)
            if (acaBlocks(row, col))
                return acaBlocks(row, col)->parallelizationOptions();
    return Fiber::ParallelizationOptions(); // actually should never get here
}

template <typename ValueType>
std::vector<unsigned int> overall_o2pDomain(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks,
        const std::vector<size_t>& colOffsets,
        size_t size)
{
    std::vector<unsigned int> o2pDomain;
    o2pDomain.reserve(size);
    for (size_t col = 0; col < acaBlocks.extent(1); ++col) {
        int chosenRow = -1;
        for (size_t row = 0; row < acaBlocks.extent(0); ++row)
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
    assert(o2pDomain.size() == size);
    return o2pDomain;
}

template <typename ValueType>
std::vector<unsigned int> overall_o2pRange(
        const Fiber::_2dArray<shared_ptr<
            const DiscreteAcaBoundaryOperator<ValueType> > >& acaBlocks,
        const std::vector<size_t>& rowOffsets,
        size_t size)
{
    std::vector<unsigned int> o2pRange;
    o2pRange.reserve(size);
    for (size_t row = 0; row < acaBlocks.extent(0); ++row) {
        int chosenCol = -1;
        for (size_t col = 0; col < acaBlocks.extent(1); ++col)
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
    assert(o2pRange.size() == size);
    return o2pRange;
}
#endif // WITH_AHMED

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

    typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
    typedef typename AcaOp::AhmedBemBlcluster AhmedBemBlcluster;
    typedef typename AcaOp::AhmedMblock AhmedMblock;
    typedef typename AcaOp::AhmedMblockArray AhmedMblockArray;
    typedef typename AcaOp::AhmedConstMblockArray AhmedConstMblockArray;

    // Convert blocks into ACA operators
    Fiber::_2dArray<shared_ptr<const AcaOp> > acaBlocks(blockRowCount,
                                                        blockColCount);
    for (size_t col = 0; col < blockColCount; ++col)
        for (size_t row = 0; row < blockRowCount; ++row)
            if (m_blocks(row, col))
                acaBlocks(row, col) = boost::shared_dynamic_cast<const AcaOp>(
                            m_blocks(row, col)->asDiscreteAcaBoundaryOperator(
                                eps, maximumRank));

    std::vector<AhmedConstMblockArray> allSharedMblocks;
    AhmedMblockArray mblocks;
    Fiber::_2dArray<size_t> indexOffsets(blockRowCount, blockColCount);
    size_t totalMblockCount = 0;
    collectMblocks(acaBlocks, m_rowCounts, m_columnCounts,
                   allSharedMblocks, mblocks, indexOffsets, totalMblockCount);

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
            for (size_t col = 0; col < blockColCount; ++col, ++i) {
                    std::auto_ptr<AhmedBemBlcluster> root(
                                new AhmedBemBlcluster(
                                    rowOffsets[row], colOffsets[col],
                                    m_rowCounts[row], m_columnCounts[col]));
                    if (acaBlocks(row, col))
                        copySonsAdjustingIndices<ValueType>(
                                    acaBlocks(row, col)->blockCluster().get(),
                                    root.get(),
                                    rowOffsets[row],
                                    colOffsets[col],
                                    indexOffsets(row, col));
                    else {
                        root->setidx(indexOffsets(row, col));
                        root->setadm(true); // unsure about it
                        root->setsep(true); // unsure about it
                    }
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
    shared_ptr<AhmedBemBlcluster> mblockCluster(
                new AhmedBemBlcluster(0, 0, totalRowCount, totalColCount));
    mblockCluster->setsons(blockRowCount, blockColCount, &roots[0]);

    // Create the domain and range permutation arrays
    std::vector<unsigned int> o2pDomain =
            overall_o2pDomain(acaBlocks, colOffsets, totalColCount);
    std::vector<unsigned int> o2pRange =
            overall_o2pRange(acaBlocks, rowOffsets, totalRowCount);

    // Gather remaining data necessary to create the combined ACA operator
    const double overallEps_ = overallEps(acaBlocks);
    const int overallMaximumRank_ = overallMaximumRank(acaBlocks);
    const int symmetry = overallSymmetry(acaBlocks);
    const Fiber::ParallelizationOptions parallelOptions =
            overallParallelizationOptions(acaBlocks);

    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new AcaOp(
                    totalRowCount, totalColCount,
                    overallEps_, overallMaximumRank_,
                    symmetry,
                    mblockCluster, mblocks,
                    o2pDomain, o2pRange,
                    parallelOptions, allSharedMblocks));
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
