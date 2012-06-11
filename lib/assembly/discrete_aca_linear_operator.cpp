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

#include "discrete_aca_linear_operator.hpp"
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
DiscreteAcaLinearOperator<ValueType>::
DiscreteAcaLinearOperator(
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
void
DiscreteAcaLinearOperator<ValueType>::
dump() const
{
    std::cout << asMatrix() << std::endl;
}

template <typename ValueType>
arma::Mat<ValueType>
DiscreteAcaLinearOperator<ValueType>::
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
DiscreteAcaLinearOperator<ValueType>::
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
DiscreteAcaLinearOperator<ValueType>::
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
DiscreteAcaLinearOperator<ValueType>::
addBlock(const std::vector<int>& rows,
         const std::vector<int>& cols,
         const ValueType alpha,
         arma::Mat<ValueType>& block) const
{
    throw std::runtime_error("DiscreteAcaLinearOperator::"
                             "addBlock(): not implemented yet");
}

template <typename ValueType>
const DiscreteAcaLinearOperator<ValueType>&
DiscreteAcaLinearOperator<ValueType>::castToAca(
        const DiscreteLinearOperator<ValueType>& discreteOperator)
{
    return dynamic_cast<const DiscreteAcaLinearOperator<ValueType>&>(discreteOperator);
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteAcaLinearOperator<ValueType>::
domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteAcaLinearOperator<ValueType>::
range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool
DiscreteAcaLinearOperator<ValueType>::
opSupportedImpl(Thyra::EOpTransp M_trans) const
{
    // TODO: implement remaining variants (transpose & conjugate transpose)
    return (M_trans == Thyra::NOTRANS);
}

template <typename ValueType>
void
DiscreteAcaLinearOperator<ValueType>::
applyImpl(const Thyra::EOpTransp M_trans,
          const Thyra::MultiVectorBase<ValueType>& X_in,
          const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> >& Y_inout,
          const ValueType alpha,
          const ValueType beta) const
{
    typedef Thyra::Ordinal Ordinal;

    // Note: the name is VERY misleading: these asserts don't disappear in
    // release runs, and in case of failure throw exceptions rather than
    // abort.
    TEUCHOS_ASSERT(X_in.range()->isCompatible(*this->domain()));
    TEUCHOS_ASSERT(Y_inout->range()->isCompatible(*this->range()));

    const Ordinal colCount = X_in.domain()->dim();

    // Loop over the input columns

    for (Ordinal col = 0; col < colCount; ++col) {
        // Get access the the elements of column col
        Thyra::ConstDetachedSpmdVectorView<ValueType> xVec(X_in.col(col));
        Thyra::DetachedSpmdVectorView<ValueType> yVec(Y_inout->col(col));
        const Teuchos::ArrayRCP<const ValueType> xArray(xVec.sv().values());
        const Teuchos::ArrayRCP<ValueType> yArray(yVec.sv().values());

        // Wrap the Trilinos array in an Armadillo vector. const_cast is used
        // because it's more natural to have a const arma::Col<ValueType> array
        // than an arma::Col<const ValueType> one.
        const arma::Col<ValueType> xCol(
                    const_cast<ValueType*>(xArray.get()), xArray.size(),
                    false /* copy_aux_mem */);
        arma::Col<ValueType> yCol(yArray.get(), yArray.size(), false);

        applyBuiltInImpl(static_cast<TranspositionMode>(M_trans),
                         xCol, yCol, alpha, beta);
    }
}
#endif // WITH_TRILINOS

template <typename ValueType>
void
DiscreteAcaLinearOperator<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
    if (trans != NO_TRANSPOSE)
        throw std::runtime_error(
                "DiscreteAcaLinearOperator::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE are not supported");
    if (columnCount() != x_in.n_rows && rowCount() != y_inout.n_rows)
        throw std::invalid_argument(
                "DiscreteAcaLinearOperator::applyBuiltInImpl(): "
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
DiscreteAcaLinearOperator<ValueType>::
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

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteAcaLinearOperator);

} // namespace Bempp

#endif // WITH_AHMED
