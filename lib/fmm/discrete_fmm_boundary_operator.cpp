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

#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_trilinos.hpp"
#ifdef WITH_AHMED

#include "discrete_fmm_boundary_operator.hpp"
#include "octree.hpp"
#include "../common/shared_ptr.hpp"

#include "../assembly/ahmed_aux.hpp"

#include "../common/chunk_statistics.hpp"
#include "../common/complex_aux.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/serial_blas_region.hpp"

#include <fstream>
#include <iostream>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/type_traits/is_complex.hpp>

#include <tbb/blocked_range.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>

#ifdef WITH_TRILINOS
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{


template <typename ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
DiscreteFmmBoundaryOperator(
		unsigned int rowCount, unsigned int columnCount,
		const shared_ptr<Octree<ValueType> > &octree,
		int symmetry_) :
#ifdef WITH_TRILINOS
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
    m_rowCount(rowCount), m_columnCount(columnCount),
#endif
	m_octree(octree),
	m_symmetry(symmetry_)
{
}

/*
template <typename ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
DiscreteFmmBoundaryOperator(
        unsigned int rowCount, unsigned int columnCount,
        double eps_,
        int maximumRank_,
        int symmetry_,
        const shared_ptr<const AhmedBemBlcluster>& blockCluster_,
        const AhmedMblockArray& blocks_,
        const IndexPermutation& domainPermutation_,
        const IndexPermutation& rangePermutation_,
        const ParallelizationOptions& parallelizationOptions_,
        const std::vector<AhmedConstMblockArray>& sharedBlocks_) :
#ifdef WITH_TRILINOS
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
    m_rowCount(rowCount), m_columnCount(columnCount),
#endif
    m_eps(eps_),
    m_maximumRank(maximumRank_),
    m_symmetry(symmetry_),
    m_blockCluster(blockCluster_), m_blocks(blocks_),
    m_domainPermutation(domainPermutation_),
    m_rangePermutation(rangePermutation_),
    m_parallelizationOptions(parallelizationOptions_),
    m_sharedBlocks(sharedBlocks_)
{
    if (eps_ <= 0 || eps_ > 1)
        std::cout << "DiscreteFmmBoundaryOperator::DiscreteFmmBoundaryOperator(): "
                     "warning: suspicious value of eps ("
                  << eps_ << ")" << std::endl;
    if (maximumRank_ < 0)
        std::cout << "DiscreteFmmBoundaryOperator::DiscreteFmmBoundaryOperator(): "
                     "warning: suspicious value of maximumRank ("
                  << maximumRank_ << ")" << std::endl;
}

template <typename ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
DiscreteFmmBoundaryOperator(
        unsigned int rowCount, unsigned int columnCount,
        int maximumRank_,
        int symmetry_,
        std::auto_ptr<const AhmedBemBlcluster> blockCluster_,
        AhmedMblockArray blocks_,
        const IndexPermutation& domainPermutation_,
        const IndexPermutation& rangePermutation_,
        const ParallelizationOptions& parallelizationOptions_,
        const std::vector<AhmedConstMblockArray>& sharedBlocks_) :
#ifdef WITH_TRILINOS
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
    m_rowCount(rowCount), m_columnCount(columnCount),
#endif
    m_eps(1e-4), // just a guess...
    m_maximumRank(maximumRank_),
    m_symmetry(symmetry_),
    m_blockCluster(blockCluster_.release()), m_blocks(blocks_),
    m_domainPermutation(domainPermutation_),
    m_rangePermutation(rangePermutation_),
    m_parallelizationOptions(parallelizationOptions_),
    m_sharedBlocks(sharedBlocks_)
{
}

template <typename ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
DiscreteFmmBoundaryOperator(
        unsigned int rowCount, unsigned int columnCount,
        int maximumRank_,
        int symmetry_,
        const shared_ptr<const AhmedBemBlcluster>& blockCluster_,
        const AhmedMblockArray& blocks_,
        const IndexPermutation& domainPermutation_,
        const IndexPermutation& rangePermutation_,
        const ParallelizationOptions& parallelizationOptions_,
        const std::vector<AhmedConstMblockArray>& sharedBlocks_) :
#ifdef WITH_TRILINOS
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
    m_rowCount(rowCount), m_columnCount(columnCount),
#endif
    m_eps(1e-4), // just a guess...
    m_maximumRank(maximumRank_),
    m_symmetry(symmetry_),
    m_blockCluster(blockCluster_), m_blocks(blocks_),
    m_domainPermutation(domainPermutation_),
    m_rangePermutation(rangePermutation_),
    m_parallelizationOptions(parallelizationOptions_),
    m_sharedBlocks(sharedBlocks_)
{
}
*/
template <typename ValueType>
arma::Mat<ValueType>
DiscreteFmmBoundaryOperator<ValueType>::
asMatrix() const
{
/*    const unsigned int nRows = rowCount();
    const unsigned int nCols = columnCount();
    const blcluster* blockCluster = m_blockCluster.get();
    blcluster* nonconstBlockCluster = const_cast<blcluster*>(blockCluster);

    arma::Mat<ValueType> permutedOutput(nRows, nCols );
    permutedOutput.fill(0.);
    arma::Col<ValueType> unit(nCols );
    unit.fill(0.);

    for (unsigned int col = 0; col < nCols ; ++col)
    {
        if (col > 0)
            unit(col - 1) = 0.;
        unit(col) = 1.;
        if (m_symmetry & SYMMETRIC)
            mltaSyHVec(1., nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(unit.memptr()),
                       ahmedCast(permutedOutput.colptr(col)));
        else if (m_symmetry & HERMITIAN)
            mltaHeHVec(1., nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(unit.memptr()),
                       ahmedCast(permutedOutput.colptr(col)));
        else
            mltaGeHVec(1., nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(unit.memptr()),
                       ahmedCast(permutedOutput.colptr(col)));
    }

    arma::Mat<ValueType> output(nRows, nCols );
    for (unsigned int col = 0; col < nCols ; ++col)
        for (unsigned int row = 0; row < nRows; ++row)
            output(row, col) =
                    permutedOutput(m_rangePermutation.permuted(row),
                                   m_domainPermutation.permuted(col));

    return output;*/
}

template <typename ValueType>
unsigned int
DiscreteFmmBoundaryOperator<ValueType>::
rowCount() const
{
#ifdef WITH_TRILINOS
    return m_rangeSpace->dim();
#else
    return m_rowCount;
#endif
}

template <typename ValueType>
unsigned int
DiscreteFmmBoundaryOperator<ValueType>::
columnCount() const
{
#ifdef WITH_TRILINOS
    return m_domainSpace->dim();
#else
    return m_columnCount;
#endif
}

template <typename ValueType>
void
DiscreteFmmBoundaryOperator<ValueType>::
addBlock(const std::vector<int>& rows,
         const std::vector<int>& cols,
         const ValueType alpha,
         arma::Mat<ValueType>& block) const
{
    throw std::runtime_error("DiscreteFmmBoundaryOperator::"
                             "addBlock(): not implemented yet");
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::asDiscreteFmmBoundaryOperator(
        double eps, int maximumRank, bool interleave) const
{
    return this->shared_from_this(); // this-> needed for template name resolution.
}

template <typename ValueType>
const DiscreteFmmBoundaryOperator<ValueType>&
DiscreteFmmBoundaryOperator<ValueType>::castToFmm(
        const DiscreteBoundaryOperator<ValueType>& discreteOperator)
{
    return dynamic_cast<const DiscreteFmmBoundaryOperator<ValueType>&>(
                discreteOperator);
}

template <typename ValueType>
shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::castToFmm(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
        discreteOperator)
{
    shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > result =
        boost::dynamic_pointer_cast<const DiscreteFmmBoundaryOperator<ValueType> >(
                discreteOperator);
    if (result.get()==0 && (discreteOperator.get()!=0)) throw std::bad_cast();
    return result;
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::
domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteFmmBoundaryOperator<ValueType>::
range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool
DiscreteFmmBoundaryOperator<ValueType>::
opSupportedImpl(Thyra::EOpTransp M_trans) const
{
    // TODO: implement remaining variant (conjugate)
    return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS ||
            M_trans == Thyra::CONJTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void
DiscreteFmmBoundaryOperator<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
/*    if (trans != NO_TRANSPOSE && trans != TRANSPOSE && trans != CONJUGATE_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteFmmBoundaryOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE, TRANSPOSE and "
                "CONJUGATE_TRANSPOSE are not supported");
    bool transposed = (trans & TRANSPOSE);

    const blcluster* blockCluster = m_blockCluster.get();
    blcluster* nonconstBlockCluster = const_cast<blcluster*>(blockCluster);

    if ((!transposed && (columnCount() != x_in.n_rows ||
                         rowCount() != y_inout.n_rows)) ||
            (transposed && (rowCount() != x_in.n_rows ||
                            columnCount() != y_inout.n_rows)))
        throw std::invalid_argument(
                "DiscreteFmmBoundaryOperator::applyBuiltInImpl(): "
                "incorrect vector length");

    if (beta == static_cast<ValueType>(0.))
        y_inout.fill(static_cast<ValueType>(0.));
    else
        y_inout *= beta;

    arma::Col<ValueType> permutedArgument;
    if (!transposed)
        m_domainPermutation.permuteVector(x_in, permutedArgument);
    else
        m_rangePermutation.permuteVector(x_in, permutedArgument);

    arma::Col<ValueType> permutedResult;
    if (!transposed)
        m_rangePermutation.permuteVector(y_inout, permutedResult);
    else
        m_domainPermutation.permuteVector(y_inout, permutedResult);

    if (m_symmetry & SYMMETRIC) {
        // TODO: parallelise
        if (trans == NO_TRANSPOSE || trans == TRANSPOSE)
            mltaSyHVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(permutedArgument.memptr()),
                       ahmedCast(permutedResult.memptr()));
        else // trans == CONJUGATE_TRANSPOSE
            mltaSyHhVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(permutedArgument.memptr()),
                       ahmedCast(permutedResult.memptr()));
    }
    else if (m_symmetry & HERMITIAN) {
        // NO_TRANSPOSE and CONJUGATE_TRANSPOSE are equivalent
        if (trans == NO_TRANSPOSE || trans == CONJUGATE_TRANSPOSE)
            mltaHeHVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(permutedArgument.memptr()),
                       ahmedCast(permutedResult.memptr()));
        else { // trans == TRANSPOSE
            // alpha A^T x + beta y = (alpha^* A^H x^* + beta^* y^*)^*
            // = (alpha^* A x^* + beta^* y^*)^*
            permutedArgument = arma::conj(permutedArgument);
            permutedResult = arma::conj(permutedResult);
            ValueType alphaConj = conj(alpha);
            mltaHeHVec(ahmedCast(alphaConj), nonconstBlockCluster, m_blocks.get(),
                       ahmedCast(permutedArgument.memptr()),
                       ahmedCast(permutedResult.memptr()));
            permutedResult = arma::conj(permutedResult);
        }
    }
    else {
//        if (trans == NO_TRANSPOSE)
//            mltaGeHVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
//                       ahmedCast(permutedArgument.memptr()),
//                       ahmedCast(permutedResult.memptr()));
//        else // trans == CONJUGATE_TRANSPOSE
//            mltaGeHhVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
//                        ahmedCast(permutedArgument.memptr()),
//                        ahmedCast(permutedResult.memptr()));

        AhmedLeafClusterArray leafClusters(nonconstBlockCluster);
        leafClusters.sortAccordingToClusterSize();
        const size_t leafClusterCount = leafClusters.size();

        int maxThreadCount = 1;
        if (!m_parallelizationOptions.isOpenClEnabled()) {
            if (m_parallelizationOptions.maxThreadCount() ==
                    ParallelizationOptions::AUTO)
                maxThreadCount = tbb::task_scheduler_init::automatic;
            else
                maxThreadCount = m_parallelizationOptions.maxThreadCount();
        }
        tbb::task_scheduler_init scheduler(maxThreadCount);

        std::vector<ChunkStatistics> chunkStats(leafClusterCount);

        typedef MblockMultiplicationLoopBody<ValueType> Body;
        typename Body::LeafClusterIndexQueue leafClusterIndexQueue;
        for (size_t i = 0; i < leafClusterCount; ++i)
            leafClusterIndexQueue.push(i);

        // std::cout << "----------------------------\nperm Arg\n" << permutedArgument;
        // std::cout << "perm Res\n" << permutedResult;
        Body body(trans,
                  alpha, permutedArgument, permutedResult,
                  leafClusters, m_blocks,
                  leafClusterIndexQueue, chunkStats);
        {
            Fiber::SerialBlasRegion region;
            tbb::parallel_reduce(tbb::blocked_range<size_t>(0, leafClusterCount),
                                 body);
        }
        permutedResult = body.m_local_y;
    }
    if (!transposed)
        m_rangePermutation.unpermuteVector(permutedResult, y_inout);
    else
        m_domainPermutation.unpermuteVector(permutedResult, y_inout);
*/
    if (trans != NO_TRANSPOSE && trans != TRANSPOSE && trans != CONJUGATE_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteFmmBoundaryOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE and "
                "CONJUGATE_TRANSPOSE are not supported");

	bool transposed = (trans & TRANSPOSE);
	arma::Col<ValueType> y_in = y_inout;
	y_inout.fill(0.0);
	m_octree->matrixVectorProduct(x_in, y_inout);
	y_inout = alpha*y_inout + beta*y_in;
	//std::cout << "transposed = " << transposed << std::endl;
	//std::cout << "alpha = " << alpha << std::endl;
	//std::cout << "beta = " << beta << std::endl;
// <tt>y_inout := alpha * trans(L) * x_in + beta * y_inout</tt>,
}

template <typename ValueType>
void
DiscreteFmmBoundaryOperator<ValueType>::
makeAllMblocksDense()
{
    for (unsigned int i = 0; i < m_blockCluster->nleaves(); ++i)
        if (m_blocks[i]->isLrM())
            m_blocks[i]->convLrM_toGeM();
}

template <typename ValueType>
double
DiscreteFmmBoundaryOperator<ValueType>::eps() const
{
    return m_eps;
}

template <typename ValueType>
int
DiscreteFmmBoundaryOperator<ValueType>::maximumRank() const
{
    return m_maximumRank;
}

template <typename ValueType>
int
DiscreteFmmBoundaryOperator<ValueType>::actualMaximumRank() const
{
    // const_cast because Ahmed is not const-correct
    return Hmax_rank(const_cast<AhmedBemBlcluster*>(m_blockCluster.get()),
                     m_blocks.get());
}

template <typename ValueType>
int
DiscreteFmmBoundaryOperator<ValueType>::symmetry() const
{
    return m_symmetry;
}

template <typename ValueType>
shared_ptr<const typename DiscreteFmmBoundaryOperator<ValueType>::AhmedBemBlcluster>
DiscreteFmmBoundaryOperator<ValueType>::blockCluster() const
{
    return m_blockCluster;
}

template <typename ValueType>
typename DiscreteFmmBoundaryOperator<ValueType>::AhmedMblockArray
DiscreteFmmBoundaryOperator<ValueType>::blocks() const
{
    return m_blocks;
}

template <typename ValueType>
size_t
DiscreteFmmBoundaryOperator<ValueType>::blockCount() const
{
    return m_blockCluster->nleaves();
}

template <typename ValueType>
const IndexPermutation&
DiscreteFmmBoundaryOperator<ValueType>::domainPermutation() const
{
    //return m_domainPermutation;
}

template <typename ValueType>
const IndexPermutation&
DiscreteFmmBoundaryOperator<ValueType>::rangePermutation() const
{
    //return m_rangePermutation;
}

template <typename ValueType>
const ParallelizationOptions&
DiscreteFmmBoundaryOperator<ValueType>::parallelizationOptions() const
{
    return m_parallelizationOptions;
}

template <typename ValueType>
std::vector<typename DiscreteFmmBoundaryOperator<ValueType>::AhmedConstMblockArray >
DiscreteFmmBoundaryOperator<ValueType>::sharedBlocks() const
{
    return m_sharedBlocks;
}

// Global routines

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > fmmOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank)
{
/*    if (!op1 || !op2)
        throw std::invalid_argument("fmmOperatorSum(): "
                                    "both operands must be non-null");
    if (eps <= 0)
        throw std::invalid_argument("fmmOperatorSum(): eps must be positive");
    if (maximumRank < 1)
        maximumRank = 1;
    shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > fmmOp1 =
            DiscreteFmmBoundaryOperator<ValueType>::castToFmm(op1);
    shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > fmmOp2 =
            DiscreteFmmBoundaryOperator<ValueType>::castToFmm(op2);
    if (fmmOp1->symmetry() != fmmOp2->symmetry())
        throw std::runtime_error("fmmOperatorSum(): addition of two H-matrices "
                                 "of different symmetry is not supported yet");
    if (!areEqual(fmmOp1->blockCluster().get(), fmmOp2->blockCluster().get()))
        throw std::invalid_argument("fmmOperatorSum(): block cluster trees of "
                                    "both operands must be identical");
    if (fmmOp1->m_domainPermutation != fmmOp2->m_domainPermutation ||
            fmmOp1->m_rangePermutation != fmmOp2->m_rangePermutation)
        throw std::invalid_argument("fmmOperatorSum(): domain and range "
                                    "index permutations of "
                                    "both operands must be identical");

    typedef typename DiscreteFmmBoundaryOperator<ValueType>::AhmedBemBlcluster
            AhmedBemBlcluster;
    typedef typename DiscreteFmmBoundaryOperator<ValueType>::AhmedMblock
            AhmedMblock;

    shared_ptr<const AhmedBemBlcluster> sumBlockCluster = fmmOp1->blockCluster();
    blcluster* nonConstSumBlockCluster =
            const_cast<blcluster*>( // AHMED is not const-correct
                static_cast<const blcluster*>( // safe -- upcast
                                        sumBlockCluster.get()));
    boost::shared_array<AhmedMblock*> sumBlocks =
            allocateAhmedMblockArray<ValueType>(sumBlockCluster.get());
    copyH(nonConstSumBlockCluster, fmmOp1->m_blocks.get(), sumBlocks.get());
    addGeHGeH(nonConstSumBlockCluster, sumBlocks.get(), fmmOp2->m_blocks.get(),
              eps, maximumRank);
    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new DiscreteFmmBoundaryOperator<ValueType> (
                    fmmOp1->rowCount(), fmmOp1->columnCount(), eps,
                    maximumRank,
                    fmmOp1->m_symmetry & fmmOp2->m_symmetry,
                    sumBlockCluster, sumBlocks,
                    fmmOp1->m_domainPermutation, fmmOp2->m_rangePermutation,
                    fmmOp1->m_parallelizationOptions));
    return result;
*/
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
/*    if (!op)
        throw std::invalid_argument("scaledFmmOperator(): "
                                    "operand must not be null");

    typedef typename DiscreteFmmBoundaryOperator<ValueType>::AhmedBemBlcluster
            AhmedBemBlcluster;
    typedef typename DiscreteFmmBoundaryOperator<ValueType>::AhmedMblock
            AhmedMblock;

    typedef typename AhmedTypeTraits<ValueType>::Type AhmedValueType;
    const AhmedValueType ahmedMultiplier = ahmedCast(multiplier);

    shared_ptr<const DiscreteFmmBoundaryOperator<ValueType> > fmmOp =
            DiscreteFmmBoundaryOperator<ValueType>::castToFmm(op);
    const int symmetry = fmmOp->symmetry();
    if (imagPart(multiplier) != 0. &&
            symmetry & HERMITIAN && !(symmetry & SYMMETRIC))
        throw std::runtime_error(
                "scaledFmmOperator(): multiplication of non-symmetric Hermitian"
                "matrices by scalars with non-zero imaginary part is not "
                "supported yet");
    shared_ptr<const AhmedBemBlcluster> scaledBlockCluster = fmmOp->blockCluster();
    blcluster* nonConstScaledBlockCluster =
            const_cast<blcluster*>( // AHMED is not const-correct
                static_cast<const blcluster*>( // safe -- upcast
                                        scaledBlockCluster.get()));
    boost::shared_array<AhmedMblock*> scaledBlocks =
            allocateAhmedMblockArray<ValueType>(scaledBlockCluster.get());
    copyH(nonConstScaledBlockCluster, fmmOp->blocks().get(), scaledBlocks.get());

    const size_t blockCount = scaledBlockCluster->nleaves();
    for (size_t b = 0; b < blockCount; ++b) {
        AhmedMblock* block = scaledBlocks[b];
        typename AhmedTypeTraits<ValueType>::Type* data = block->getdata();
        if (block->isLrM()) // low-rank
            for (size_t i = 0; i < block->getn1() * block->rank(); ++i)
                data[i] *= ahmedMultiplier;
        else {
             // we don't support complex Hermitian blocks yet
            assert(!(block->isHeM() && boost::is_complex<ValueType>()));
            if (!block->isLtM()) // dense, but not lower-triangular
                for (size_t i = 0; i < block->nvals(); ++i)
                    data[i] *= ahmedMultiplier;
            else { // lower-triangular
                // diagonal entries have special meaning and should not be rescaled
                size_t size = block->getn1();
                assert(block->getn2() == size);
                for (size_t c = 0; c < size; ++c)
                    for (size_t r = 1; r < size; ++r)
                        data[c * size + r] *= ahmedMultiplier;
            }
        }
    }

    int scaledSymmetry = symmetry;
    if (imagPart(multiplier) != 0.)
        scaledSymmetry &= ~HERMITIAN;
    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new DiscreteFmmBoundaryOperator<ValueType> (
                    fmmOp->rowCount(), fmmOp->columnCount(),
                    fmmOp->eps(), fmmOp->maximumRank(),
                    scaledSymmetry,
                    scaledBlockCluster, scaledBlocks,
                    fmmOp->domainPermutation(), fmmOp->rangePermutation(),
                    fmmOp->parallelizationOptions()));
    return result;
*/
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledFmmOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier)
{
    return scaledFmmOperator(multiplier, op);
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteFmmBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(RESULT) \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        fmmOperatorSum( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op1, \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op2, \
            double eps, int maximumRank); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        scaledFmmOperator( \
            const RESULT& multiplier, \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        scaledFmmOperator( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op, \
            const RESULT& multiplier)

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && (defined(ENABLE_COMPLEX_BASIS_FUNCTIONS) || defined(ENABLE_COMPLEX_KERNELS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<float>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && (defined(ENABLE_COMPLEX_BASIS_FUNCTIONS) || defined(ENABLE_COMPLEX_KERNELS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>);
#endif

} // namespace Bempp

#endif // WITH_AHMED
