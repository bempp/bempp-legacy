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

#include "discrete_aca_boundary_operator.hpp"
#include "../common/shared_ptr.hpp"

#include "ahmed_aux.hpp"
#include "aca_approximate_lu_inverse.hpp"

#include "../common/complex_aux.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <fstream>
#include <iostream>
#include <boost/smart_ptr/shared_ptr.hpp>

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

namespace
{

bool areEqual(const blcluster* op1, const blcluster* op2)
{
    if (!op1 || !op2)
        return (!op1 && !op2);
    if (op1->getnrs() != op2->getnrs() || op1->getncs() != op2->getncs())
        return false;
    if (op1->isleaf())
        return (op1->getb1() == op2->getb1() && op1->getb2() == op2->getb2() &&
                op1->getn1() == op2->getn1() && op1->getn2() == op2->getn2());
    for (unsigned int rs = 0; rs < op1->getnrs(); ++rs)
        for (unsigned int cs = 0; cs < op1->getncs(); ++cs)
            if (!areEqual(op1->getson(rs, cs), op2->getson(rs, cs)))
                return false;
    return true;
}

} // namespace

template <typename ValueType>
DiscreteAcaBoundaryOperator<ValueType>::
DiscreteAcaBoundaryOperator(
        unsigned int rowCount, unsigned int columnCount,
        int maximumRank,
        Symmetry symmetry,
        std::auto_ptr<AhmedBemBlcluster> blockCluster,
        boost::shared_array<AhmedMblock*> blocks,
        const IndexPermutation& domainPermutation,
        const IndexPermutation& rangePermutation) :
#ifdef WITH_TRILINOS
    m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
    m_rowCount(rowCount), m_columnCount(columnCount),
#endif
    m_maximumRank(maximumRank),
    m_symmetry(symmetry),
    m_blockCluster(blockCluster), m_blocks(blocks),
    m_domainPermutation(domainPermutation),
    m_rangePermutation(rangePermutation)
{
}

template <typename ValueType>
arma::Mat<ValueType>
DiscreteAcaBoundaryOperator<ValueType>::
asMatrix() const
{
    const unsigned int nRows = rowCount();
    const unsigned int nCols = columnCount();

    arma::Mat<ValueType> permutedOutput(nRows, nCols );
    permutedOutput.fill(0.);
    arma::Col<ValueType> unit(nCols );
    unit.fill(0.);

    for (unsigned int col = 0; col < nCols ; ++col)
    {
        if (col > 0)
            unit(col - 1) = 0.;
        unit(col) = 1.;
        if (m_symmetry & HERMITIAN)
#ifdef AHMED_PRERELEASE
            multaHSymvec
#else
            mltaHeHVec
#endif
                (1., m_blockCluster.get(), m_blocks.get(),
                      ahmedCast(unit.memptr()),
                      ahmedCast(permutedOutput.colptr(col)));
        else
#ifdef AHMED_PRERELEASE
            multaHvec
#else
            mltaGeHVec
#endif
                (1., m_blockCluster.get(), m_blocks.get(),
                      ahmedCast(unit.memptr()),
                      ahmedCast(permutedOutput.colptr(col)));
    }

    arma::Mat<ValueType> output(nRows, nCols );
    for (unsigned int col = 0; col < nCols ; ++col)
        for (unsigned int row = 0; row < nRows; ++row)
            output(row, col) =
                    permutedOutput(m_rangePermutation.permuted(row),
                                   m_domainPermutation.permuted(col));

    return output;
}

template <typename ValueType>
unsigned int
DiscreteAcaBoundaryOperator<ValueType>::
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
DiscreteAcaBoundaryOperator<ValueType>::
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
DiscreteAcaBoundaryOperator<ValueType>::
addBlock(const std::vector<int>& rows,
         const std::vector<int>& cols,
         const ValueType alpha,
         arma::Mat<ValueType>& block) const
{
    throw std::runtime_error("DiscreteAcaBoundaryOperator::"
                             "addBlock(): not implemented yet");
}

template <typename ValueType>
const DiscreteAcaBoundaryOperator<ValueType>&
DiscreteAcaBoundaryOperator<ValueType>::castToAca(
        const DiscreteBoundaryOperator<ValueType>& discreteOperator)
{
    return dynamic_cast<const DiscreteAcaBoundaryOperator<ValueType>&>(
                discreteOperator);
}

template <typename ValueType>
shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> >
DiscreteAcaBoundaryOperator<ValueType>::castToAca(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >&
        discreteOperator)
{
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > result =
        boost::dynamic_pointer_cast<const DiscreteAcaBoundaryOperator<ValueType> >(
                discreteOperator);
    if (result.get()==0 && (discreteOperator.get()!=0)) throw std::bad_cast();
    return result;
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteAcaBoundaryOperator<ValueType>::
domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteAcaBoundaryOperator<ValueType>::
range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool
DiscreteAcaBoundaryOperator<ValueType>::
opSupportedImpl(Thyra::EOpTransp M_trans) const
{
    // TODO: implement remaining variants (transpose & conjugate)
    return (M_trans == Thyra::NOTRANS || M_trans == Thyra::CONJTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void
DiscreteAcaBoundaryOperator<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
#ifdef AHMED_PRERELEASE
    if (trans != NO_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE are not supported");
#else
    if (trans != NO_TRANSPOSE && trans != CONJUGATE_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE and "
                "CONJUGATE_TRANSPOSE are not supported");
#endif
    bool transposed = trans == CONJUGATE_TRANSPOSE;

    if ((!transposed && (columnCount() != x_in.n_rows ||
                         rowCount() != y_inout.n_rows)) ||
            (transposed && (rowCount() != x_in.n_rows ||
                            columnCount() != y_inout.n_rows)))
        throw std::invalid_argument(
                "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
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

    if (m_symmetry & HERMITIAN) {
#ifdef AHMED_PRERELEASE
        multaHSymvec(ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                     ahmedCast(permutedArgument.memptr()),
                     ahmedCast(permutedResult.memptr()));
#else
        // NO_TRANSPOSE and CONJUGATE_TRANSPOSE are equivalent
        mltaHeHVec(ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                   ahmedCast(permutedArgument.memptr()),
                   ahmedCast(permutedResult.memptr()));
#endif
    }
    else {
#ifdef AHMED_PRERELEASE
        multaHvec(ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                  ahmedCast(permutedArgument.memptr()),
                  ahmedCast(permutedResult.memptr()));
#else
        if (trans == NO_TRANSPOSE)
            mltaGeHVec(ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                       ahmedCast(permutedArgument.memptr()),
                       ahmedCast(permutedResult.memptr()));
        else // trans == CONJUGATE_TRANSPOSE
            mltaGeHhVec(ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                        ahmedCast(permutedArgument.memptr()),
                        ahmedCast(permutedResult.memptr()));
#endif
    }
    if (!transposed)
        m_rangePermutation.unpermuteVector(permutedResult, y_inout);
    else
        m_domainPermutation.unpermuteVector(permutedResult, y_inout);
}

template <typename ValueType>
void
DiscreteAcaBoundaryOperator<ValueType>::
makeAllMblocksDense()
{
    for (unsigned int i = 0; i < m_blockCluster->nleaves(); ++i)
#ifdef AHMED_PRERELEASE
        if (m_blocks[i]->islwr())
            m_blocks[i]->conv_lwr_to_dns();
#else
        if (m_blocks[i]->isLrM())
            m_blocks[i]->convLrM_toGeM();
#endif
}

// Global routines

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorSum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2,
        double eps, int maximumRank)
{
    if (!op1 || !op2)
        throw std::invalid_argument("acaOperatorSum(): "
                                    "both operands must be non-null");
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > acaOp1 =
            DiscreteAcaBoundaryOperator<ValueType>::castToAca(op1);
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > acaOp2 =
            DiscreteAcaBoundaryOperator<ValueType>::castToAca(op2);
    if (!areEqual(acaOp1->m_blockCluster.get(), acaOp2->m_blockCluster.get()))
        throw std::invalid_argument("acaOperatorSum(): block cluster trees of "
                                    "both operands must be identical");
    if (acaOp1->m_domainPermutation != acaOp2->m_domainPermutation ||
            acaOp1->m_rangePermutation != acaOp2->m_rangePermutation)
        throw std::invalid_argument("acaOperatorSum(): domain and range "
                                    "index permutations of "
                                    "both operands must be identical");

    typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster
            AhmedBemBlcluster;
    typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblock
            AhmedMblock;

    std::auto_ptr<AhmedBemBlcluster> sumBlockCluster(
                new AhmedBemBlcluster(acaOp1->m_blockCluster.get()));
    boost::shared_array<AhmedMblock*> sumBlocks =
            allocateAhmedMblockArray<ValueType>(sumBlockCluster.get());
    copyH(sumBlockCluster.get(), acaOp1->m_blocks.get(), sumBlocks.get());
    addGeHGeH(sumBlockCluster.get(), sumBlocks.get(), acaOp2->m_blocks.get(),
              eps, maximumRank);
    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new DiscreteAcaBoundaryOperator<ValueType> (
                    acaOp1->rowCount(), acaOp1->columnCount(), maximumRank,
                    Symmetry(acaOp1->m_symmetry & acaOp2->m_symmetry),
                    sumBlockCluster, sumBlocks,
                    acaOp1->m_domainPermutation, acaOp2->m_rangePermutation));
    return result;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator(
        const ValueType& multiplier,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    if (!op)
        throw std::invalid_argument("scaledAcaOperator(): "
                                    "operand must not be null");

    typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster
            AhmedBemBlcluster;
    typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblock
            AhmedMblock;

    typedef typename AhmedTypeTraits<ValueType>::Type AhmedValueType;
    const AhmedValueType ahmedMultiplier = ahmedCast(multiplier);

    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > acaOp =
            DiscreteAcaBoundaryOperator<ValueType>::castToAca(op);
    std::auto_ptr<AhmedBemBlcluster> scaledBlockCluster(
                new AhmedBemBlcluster(acaOp->m_blockCluster.get()));
    boost::shared_array<AhmedMblock*> scaledBlocks =
            allocateAhmedMblockArray<ValueType>(scaledBlockCluster.get());
    copyH(scaledBlockCluster.get(), acaOp->m_blocks.get(), scaledBlocks.get());

    const size_t blockCount = scaledBlockCluster->nleaves();
    for (size_t b = 0; b < blockCount; ++b) {
        AhmedMblock* block = scaledBlocks[b];
        typename AhmedTypeTraits<ValueType>::Type* data = block->getdata();
        if (block->isLrM()) // low-rank
            for (size_t i = 0; i < block->getn1() * block->rank(); ++i)
                data[i] *= ahmedMultiplier;
        else if (!block->isLtM()) // dense, but not lower-triangular
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

    unsigned int scaledSymmetry = acaOp->m_symmetry;
    if (imagPart(multiplier) != 0.)
        scaledSymmetry &= ~HERMITIAN;
    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new DiscreteAcaBoundaryOperator<ValueType> (
                    acaOp->rowCount(), acaOp->columnCount(), acaOp->m_maximumRank,
                    Symmetry(scaledSymmetry),
                    scaledBlockCluster, scaledBlocks,
                    acaOp->m_domainPermutation, acaOp->m_rangePermutation));
    return result;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > scaledAcaOperator(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        const ValueType& multiplier)
{
    return scaledAcaOperator(multiplier, op);
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > acaOperatorApproximateLuInverse(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        double delta)
{
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType> > acaOp =
            DiscreteAcaBoundaryOperator<ValueType>::castToAca(op);
    shared_ptr<const DiscreteBoundaryOperator<ValueType> > result(
                new AcaApproximateLuInverse<ValueType>(*acaOp, delta));
    return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteAcaBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(RESULT) \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        acaOperatorSum( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op1, \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op2, \
            double eps, int maximumRank); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        scaledAcaOperator( \
            const RESULT& multiplier, \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        scaledAcaOperator( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op, \
            const RESULT& multiplier); \
    template shared_ptr<const DiscreteBoundaryOperator<RESULT> > \
        acaOperatorApproximateLuInverse( \
            const shared_ptr<const DiscreteBoundaryOperator<RESULT> >& op, \
            double delta)

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
