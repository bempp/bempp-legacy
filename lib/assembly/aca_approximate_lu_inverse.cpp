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
#include "aca_approximate_lu_inverse.hpp"

#include "ahmed_aux.hpp"
#include "discrete_aca_boundary_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

template <typename ValueType>
AcaApproximateLuInverse<ValueType>::AcaApproximateLuInverse(
        const DiscreteAcaBoundaryOperator<ValueType>& fwdOp,
        MagnitudeType delta) :
    // All range-domain swaps intended!
#ifdef WITH_TRILINOS
    m_domainSpace(fwdOp.m_rangeSpace),
    m_rangeSpace(fwdOp.m_domainSpace),
#else
    m_rowCount(fwdOp.columnCount()), m_columnCount(fwdOp.rowCount()),
#endif
    m_blockCluster(0), m_blocksL(0), m_blocksU(0),
    m_domainPermutation(fwdOp.m_rangePermutation),
    m_rangePermutation(fwdOp.m_domainPermutation)
{
    const blcluster* fwdBlockCluster = fwdOp.m_blockCluster.get();
    bool result = genLUprecond(const_cast<blcluster*>(fwdBlockCluster),
                               fwdOp.m_blocks.get(),
                               delta, fwdOp.m_maximumRank,
                               m_blockCluster, m_blocksL, m_blocksU, true);
    if (!result)
        throw std::runtime_error(
                "AcaApproximateLuInverse::AcaApproximateLuInverse(): "
                "Approximate LU factorisation failed");
}

template <>
AcaApproximateLuInverse<std::complex<float> >::AcaApproximateLuInverse(
        const DiscreteAcaBoundaryOperator<std::complex<float> >& fwdOp,
        MagnitudeType delta) :
    // All range-domain swaps intended!
#ifdef WITH_TRILINOS
    m_domainSpace(fwdOp.m_rangeSpace),
    m_rangeSpace(fwdOp.m_domainSpace),
#else
    m_rowCount(fwdOp.columnCount()), m_columnCount(fwdOp.rowCount()),
#endif
    m_blockCluster(0), m_blocksL(0), m_blocksU(0),
    m_domainPermutation(fwdOp.m_rangePermutation),
    m_rangePermutation(fwdOp.m_domainPermutation)
{
    // Ahmed doesn't define the genLUprecond() variant with
    // the second parameter of type mblock<scomp>**
    throw std::runtime_error(
                "AcaApproximateLuInverse::AcaApproximateLuInverse(): "
                "due to a deficiency in Ahmed approximate LU factorisation "
                "of single-precision complex H matrices is not supported");
}

template <typename ValueType>
AcaApproximateLuInverse<ValueType>::~AcaApproximateLuInverse()
{
    if (m_blockCluster)
    {
        freembls(m_blockCluster, m_blocksL);
        freembls(m_blockCluster, m_blocksU);
        delete m_blockCluster;
    }
}

template <typename ValueType>
unsigned int AcaApproximateLuInverse<ValueType>::rowCount() const
{
#ifdef WITH_TRILINOS
    return m_rangeSpace->dim();
#else
    return m_rowCount;
#endif
}

template <typename ValueType>
unsigned int AcaApproximateLuInverse<ValueType>::columnCount() const
{
#ifdef WITH_TRILINOS
    return m_domainSpace->dim();
#else
    return m_columnCount;
#endif
}

template <typename ValueType>
void AcaApproximateLuInverse<ValueType>::addBlock(
        const std::vector<int>& rows, const std::vector<int>& cols, const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    throw std::runtime_error("AcaApproximateLuInverse::addBlock(): "
                             "not implemented");
}

#ifdef WITH_TRILINOS

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
AcaApproximateLuInverse<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
AcaApproximateLuInverse<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool AcaApproximateLuInverse<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    // TODO: implement remaining variants (transpose & conjugate transpose)
    return (M_trans == Thyra::NOTRANS);
}

template <typename ValueType>
void AcaApproximateLuInverse<ValueType>::applyImpl(
        const Thyra::EOpTransp M_trans,
        const Thyra::MultiVectorBase<ValueType>& X_in,
        const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> >& Y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    typedef Thyra::Ordinal Ordinal;

    // Note: the name is VERY misleading: these asserts don't disappear in
    // release runs, and in case of failure throw exceptions rather than
    // abort.
    TEUCHOS_ASSERT(this->opSupported(M_trans));
    TEUCHOS_ASSERT(X_in.range()->isCompatible(*this->domain()));
    TEUCHOS_ASSERT(Y_inout->range()->isCompatible(*this->range()));
    TEUCHOS_ASSERT(Y_inout->domain()->isCompatible(*X_in.domain()));

    const Ordinal colCount = X_in.domain()->dim();

    // Loop over the input columns

    for (Ordinal col = 0; col < colCount; ++col) {
        // Get access the the elements of column col
        Thyra::ConstDetachedSpmdVectorView<ValueType> xVec(X_in.col(col));
        Thyra::DetachedSpmdVectorView<ValueType> yVec(Y_inout->col(col));
        const Teuchos::ArrayRCP<const ValueType> xArray(xVec.sv().values());
        const Teuchos::ArrayRCP<ValueType> yArray(yVec.sv().values());

        // const_cast because it's more natural to have
        // a const arma::Col<ValueType> array than
        // an arma::Col<const ValueType> one.
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
void AcaApproximateLuInverse<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
    if (trans != NO_TRANSPOSE)
        throw std::runtime_error(
                "AcaApproximateLuInverse::applyBuiltInImpl(): "
                "transposition modes other than NO_TRANSPOSE are not supported");
    if (columnCount() != x_in.n_rows && rowCount() != y_inout.n_rows)
        throw std::invalid_argument(
                "AcaApproximateLuInverse::applyBuiltInImpl(): "
                "incorrect vector length");

    if (beta == static_cast<ValueType>(0.))
        y_inout.fill(static_cast<ValueType>(0.));
    else
        y_inout *= beta;

    // will act both as a permuted argument and permuted result
    arma::Col<ValueType> permuted;
    m_domainPermutation.permuteVector(x_in, permuted);

    HLU_solve(m_blockCluster, m_blocksL, m_blocksU, ahmedCast(permuted.memptr()));

    arma::Col<ValueType> operatorActionResult;
    m_rangePermutation.unpermuteVector(permuted, operatorActionResult);
    y_inout += alpha * operatorActionResult;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(AcaApproximateLuInverse);

} // namespace Bempp

#endif // WITH_AHMED
