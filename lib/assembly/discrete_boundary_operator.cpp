#include "discrete_boundary_operator.hpp"
#include "discrete_boundary_operator_sum.hpp"
#include "scaled_discrete_boundary_operator.hpp"
#include "../common/shared_ptr.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include <Thyra_DetachedSpmdVectorView.hpp>

namespace Bempp
{

template <typename ValueType>
arma::Mat<ValueType> DiscreteBoundaryOperator<ValueType>::asMatrix() const
{
    // Default brute-force implementation: apply operator to all basis vectors
    const size_t nRows = rowCount();
    const size_t nCols = columnCount();
    arma::Col<ValueType> unit(nCols);
    arma::Mat<ValueType> result(nRows, nCols);
    result.fill(0.); // for safety, in case there was a bug in the handling of
                     // beta == 0. in a particular subclass' applyBuiltInImpl()
                     // override...
    unit.fill(0.);
    for (size_t i = 0; i < nCols; ++i) {
        arma::Col<ValueType> activeCol(result.unsafe_col(i));
        if (i > 0)
            unit(i - 1) = 0.;
        unit(i) = 1.;        
        applyBuiltInImpl(NO_TRANSPOSE, unit, activeCol, 1., 0.);
    }

    return result;
}

template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::dump() const
{
    std::cout << asMatrix() << std::endl;
}

#ifdef WITH_TRILINOS
template <typename ValueType>
void DiscreteBoundaryOperator<ValueType>::applyImpl(
        const Thyra::EOpTransp M_trans,
        const Thyra::MultiVectorBase<ValueType> &X_in,
        const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> > &Y_inout,
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
        // Get access the the elements of X_in's and Y_inout's column #col
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
#endif

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator+(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return op;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator-(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return static_cast<ValueType>(-1.) * op;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator+(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<const DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorSum<ValueType>(op1,op2));
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator-(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<const DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorSum<ValueType>(op1,static_cast<ValueType>(-1.)*op2));
}


template <typename ValueType, typename ScalarType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator*(
        ScalarType scalar,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return  shared_ptr<const DiscreteBoundaryOperator<ValueType> >(
                new ScaledDiscreteBoundaryOperator<ValueType>(static_cast<ValueType>(scalar),op));
}

template <typename ValueType, typename ScalarType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator*(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        ScalarType scalar)
{
    return shared_ptr<const DiscreteBoundaryOperator<ValueType> >(
                new ScaledDiscreteBoundaryOperator<ValueType>(static_cast<ValueType>(scalar),op));
}

template <typename ValueType, typename ScalarType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator/(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        ScalarType scalar)
{
    if (scalar == static_cast<ScalarType>(0.))
        throw std::runtime_error("operator/(DiscreteBoundaryOperator, scalar): "
                                     "Division by zero");
    return shared_ptr<const DiscreteBoundaryOperator<ValueType> >(
                new ScaledDiscreteBoundaryOperator<ValueType>(
                    static_cast<ValueType>(static_cast<ScalarType>(1.)/scalar),op));
}


FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(VALUE) \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator+( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator-( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator+( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2); \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator-( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2);
#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR( VALUE , SCALAR ) \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator*( \
        SCALAR scalar, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator*( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op,\
        SCALAR scalar); \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator/( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op, \
        SCALAR scalar);

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, float );
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, double);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && (defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && (defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, std::complex<double>);
#endif



} // namespace Bempp
