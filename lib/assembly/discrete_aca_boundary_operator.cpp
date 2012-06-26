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

#include "config_ahmed.hpp"
#include "config_trilinos.hpp"
#ifdef WITH_AHMED

#include "discrete_aca_boundary_operator.hpp"
#include "ahmed_aux.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <fstream>

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

template <typename ValueType>
DiscreteAcaBoundaryOperator<ValueType>::
DiscreteAcaBoundaryOperator(
        unsigned int rowCount, unsigned int columnCount,
        int maximumRank,
        bool symmetric,
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
    m_symmetric(symmetric),
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
        if (m_symmetric)
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
    return dynamic_cast<const DiscreteAcaBoundaryOperator<ValueType>&>(discreteOperator);
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
    // TODO: implement remaining variants (transpose & conjugate transpose)
    return (M_trans == Thyra::NOTRANS);
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
    if (trans != NO_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE are not supported");
    if (columnCount() != x_in.n_rows || rowCount() != y_inout.n_rows)
        throw std::invalid_argument(
                "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
                "incorrect vector length");

    if (beta == static_cast<ValueType>(0.))
        y_inout.fill(static_cast<ValueType>(0.));
    else
        y_inout *= beta;

    arma::Col<ValueType> permutedArgument;
    m_domainPermutation.permuteVector(x_in, permutedArgument);

    // const_cast because Ahmed internally calls BLAS
    // functions, which don't respect const-correctness
    arma::Col<ValueType> permutedResult;
    m_rangePermutation.permuteVector(y_inout, permutedResult);
    if (m_symmetric)
#ifdef AHMED_PRERELEASE
            multaHSymvec
#else
            mltaHeHVec
#endif
                (ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                     ahmedCast(permutedArgument.memptr()),
                     ahmedCast(permutedResult.memptr()));
    else
#ifdef AHMED_PRERELEASE
            multaHvec
#else
            mltaGeHVec
#endif
                (ahmedCast(alpha), m_blockCluster.get(), m_blocks.get(),
                  ahmedCast(permutedArgument.memptr()),
                  ahmedCast(permutedResult.memptr()));
    m_rangePermutation.unpermuteVector(permutedResult, y_inout);
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

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteAcaBoundaryOperator);

} // namespace Bempp

#endif // WITH_AHMED
