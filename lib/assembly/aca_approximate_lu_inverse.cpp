#include "../common/config_ahmed.hpp"
#include "../common/config_trilinos.hpp"

#ifdef WITH_AHMED
#include "aca_approximate_lu_inverse.hpp"

#include "ahmed_aux.hpp"
#include "discrete_aca_scalar_valued_linear_operator.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

template <typename ValueType>
AcaApproximateLuInverse<ValueType>::AcaApproximateLuInverse(
        const DiscreteAcaScalarValuedLinearOperator<ValueType>& fwdOp,
        ValueType delta) :
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
    bool result = genLUprecond(fwdOp.m_blockCluster.get(), fwdOp.m_blocks.get(),
                               delta, fwdOp.m_maximumRank,
                               m_blockCluster, m_blocksL, m_blocksU, true);
    if (!result)
        throw std::runtime_error(
                "AcaApproximateLuInverse::AcaApproximateLuInverse(): "
                "Approximate LU factorisation failed");
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
void AcaApproximateLuInverse<ValueType>::dump() const
{
    throw std::runtime_error("AcaApproximateLuInverse::dump(): "
                             "not implemented");
}

template <typename ValueType>
arma::Mat<ValueType> AcaApproximateLuInverse<ValueType>::asMatrix() const
{
    // doesn't make sense for this operator
    throw std::runtime_error("AcaApproximateLuInverse::asMatrix(): "
                             "not implemented");
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
        const std::vector<int>& rows, const std::vector<int>& cols,
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
    //    std::cout << "alpha: " << alpha << "\n";
    //    std::cout << "beta: " << beta << "\n";

    typedef Thyra::Ordinal Ordinal;

    // Note: the name is VERY misleading: these asserts don't disappear in
    // release runs, and in case of failure throw exceptions rather than
    // abort.
    TEUCHOS_ASSERT(this->opSupported(M_trans));
    TEUCHOS_ASSERT(X_in.range()->isCompatible(*this->domain()));
    TEUCHOS_ASSERT(Y_inout->range()->isCompatible(*this->range()));
    TEUCHOS_ASSERT(Y_inout->domain()->isCompatible(*X_in.domain()));

    //    std::cout << "Matrix top-left:\n" << asMatrix()(arma::span(0,3), arma::span(0,3));

    const Ordinal colCount = X_in.domain()->dim();

    // Loop over the input columns

    for (Ordinal col = 0; col < colCount; ++col) {
        // Get access the the elements of column col
        Thyra::ConstDetachedSpmdVectorView<ValueType> xVec(X_in.col(col));
        Thyra::DetachedSpmdVectorView<ValueType> yVec(Y_inout->col(col));
        const Teuchos::ArrayRCP<const ValueType> xArray(xVec.sv().values());
        const Teuchos::ArrayRCP<ValueType> yArray(yVec.sv().values());

        //        std::cout << "\n\nxArray:\n" << xArray << std::endl;
        //        for (int i = xArray.lowerOffset(); i <= xArray.upperOffset(); ++i)
        //            std::cout << xArray[i] << ", ";
        //        std::cout << "\n";
        //        std::cout << "\nyArray:\n" << yArray << std::endl;
        //        for (int i = yArray.lowerOffset(); i <= yArray.upperOffset(); ++i)
        //            std::cout << yArray[i] << ", ";
        //        std::cout << "\n";

        // const_cast because it's more natural to have
        // a const arma::Col<ValueType> array than
        // an arma::Col<const ValueType> one.
        const arma::Col<ValueType> xCol(
                    const_cast<ValueType*>(xArray.get()), xArray.size(),
                    false /* copy_aux_mem */);
        arma::Col<ValueType> yCol(yArray.get(), yArray.size(), false);

        //        std::cout << "xCol:\n" << xCol.t();
        //        std::cout << "yCol:\n" << yCol.t();

        applyBuiltInImpl(static_cast<TranspositionMode>(M_trans),
                         xCol, yCol, alpha, beta);

        //        std::cout << "\nmodified yArray:\n" << yArray << std::endl;
        //        for (int i = yArray.lowerOffset(); i <= yArray.upperOffset(); ++i)
        //            std::cout << yArray[i] << ", ";
        //        std::cout << "\n";
        //        std::cout << "modified yCol:\n" << yCol.t();

        //        if (alpha == 1. && beta == 0.)
        //            std::cout << "\n\n\nDIAGNOSTICS\n\n"
        //                      << "xCol\n" << xCol << "\n"
        //                      << "mat\n" << asMatrix() << "\n"
        //                      << "yCol\n" << yCol << "\n\n\n\n########\n";

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

    if (beta == 0.)
        y_inout.fill(0.);
    else
        y_inout *= beta;

    // will act both as a permuted argument and permuted result
    arma::Col<ValueType> permuted;
    m_domainPermutation.permuteVector(x_in, permuted);

    HLU_solve(m_blockCluster, m_blocksL, m_blocksU, permuted.memptr());

    arma::Col<ValueType> operatorActionResult;
    m_rangePermutation.unpermuteVector(permuted, operatorActionResult);
    y_inout += alpha * operatorActionResult;
}



#ifdef COMPILE_FOR_FLOAT
template class AcaApproximateLuInverse<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class AcaApproximateLuInverse<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class AcaApproximateLuInverse<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class AcaApproximateLuInverse<std::complex<double> >;
#endif

} // namespace Bempp

#endif // WITH_AHMED
