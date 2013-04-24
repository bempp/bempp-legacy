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

#include "discrete_boundary_operator.hpp"

#include "complexified_discrete_boundary_operator.hpp"
#include "discrete_boundary_operator_sum.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "scaled_discrete_boundary_operator.hpp"
#include "transposed_discrete_boundary_operator.hpp"
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
void
DiscreteBoundaryOperator<ValueType>::apply(
        const TranspositionMode trans,
        const arma::Mat<ValueType>& x_in,
        arma::Mat<ValueType>& y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    bool transposed = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);
    if (x_in.n_rows != (transposed ? rowCount() : columnCount()))
        throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                    "vector x_in has invalid length");
    if (y_inout.n_rows != (transposed ? columnCount() : rowCount()))
        throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                    "vector y_inout has invalid length");
    if (x_in.n_cols != y_inout.n_cols)
        throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                    "vectors x_in and y_inout must have "
                                    "the same number of columns");

    for (size_t i = 0; i < x_in.n_cols; ++i) {
        const arma::Col<ValueType> x_in_col = x_in.unsafe_col(i);
        arma::Col<ValueType> y_inout_col = y_inout.unsafe_col(i);
        applyBuiltInImpl(trans, x_in_col, y_inout_col, alpha, beta);
    }
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> >
DiscreteBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
        double eps, int maximumRank, bool interleave) const
{
#ifdef WITH_AHMED
    throw std::runtime_error("DiscreteBoundaryOperator::"
                             "asDiscreteAcaBoundaryOperator(): "
                             "not implemented for operators of class " +
                             std::string(typeid(*this).name()) +
                             ".");
#else
    throw std::runtime_error("DiscreteBoundaryOperator::"
                             "asDiscreteAcaBoundaryOperator(): "
                             "ACA operators are not supported because BEM++ "
                             "has been compiled without AHMED.");
#endif
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
shared_ptr<DiscreteBoundaryOperator<ValueType> > operator+(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorSum<ValueType>(op1, op2));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > sum(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorSum<ValueType>(op1, op2));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > operator-(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorSum<ValueType>(
                    op1, static_cast<ValueType>(-1.) * op2));
}

template <typename ValueType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    shared_ptr<DiscreteBoundaryOperator<ValueType> > >::type
operator*(
        ScalarType scalar,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new ScaledDiscreteBoundaryOperator<ValueType>(
                    static_cast<ValueType>(scalar), op));
}

template <typename ValueType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    shared_ptr<DiscreteBoundaryOperator<ValueType> > >::type
operator*(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        ScalarType scalar)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new ScaledDiscreteBoundaryOperator<ValueType>(
                    static_cast<ValueType>(scalar), op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > operator*(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorComposition<ValueType>(
                    op1, op2));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > mul(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new DiscreteBoundaryOperatorComposition<ValueType>(
                    op1, op2));
}

template <typename ValueType, typename ScalarType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > operator/(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        ScalarType scalar)
{
    if (scalar == static_cast<ScalarType>(0.))
        throw std::runtime_error("operator/(DiscreteBoundaryOperator, scalar): "
                                     "Division by zero");
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
                new ScaledDiscreteBoundaryOperator<ValueType>(
                    static_cast<ValueType>(static_cast<ScalarType>(1.)/scalar),op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > transpose(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
        new TransposedDiscreteBoundaryOperator<ValueType>(
            TRANSPOSE, op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > conjugate(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
        new TransposedDiscreteBoundaryOperator<ValueType>(
            CONJUGATE, op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > conjugateTranspose(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
        new TransposedDiscreteBoundaryOperator<ValueType>(
            CONJUGATE_TRANSPOSE, op));
}

template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > transpose(
        TranspositionMode trans,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op)
{
    return shared_ptr<DiscreteBoundaryOperator<ValueType> >(
        new TransposedDiscreteBoundaryOperator<ValueType>(
            trans, op));
}

template <typename RealType>
shared_ptr<DiscreteBoundaryOperator<std::complex<RealType> > > complexify(
        const shared_ptr<const DiscreteBoundaryOperator<RealType> >& op)
{
    return shared_ptr<DiscreteBoundaryOperator<std::complex<RealType> > >(
        new ComplexifiedDiscreteBoundaryOperator<RealType>(op));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(VALUE) \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator+( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<const DiscreteBoundaryOperator< VALUE > > operator-( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > operator+( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > sum( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > operator-( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > operator*( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > mul( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op1, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op2); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > transpose( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > conjugate( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > conjugateTranspose( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > transpose( \
        TranspositionMode trans, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op);

#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR( VALUE , SCALAR ) \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > operator*( \
        SCALAR scalar, \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > operator*( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op,\
        SCALAR scalar); \
    template shared_ptr<DiscreteBoundaryOperator< VALUE > > operator/( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op, \
        SCALAR scalar);

#define INSTANTIATE_FREE_FUNCTIONS_REAL_ONLY(VALUE) \
    template shared_ptr<DiscreteBoundaryOperator< std::complex<VALUE> > > complexify( \
        const shared_ptr<const DiscreteBoundaryOperator< VALUE > >& op);

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, float );
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, double);
INSTANTIATE_FREE_FUNCTIONS_REAL_ONLY(float);
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
INSTANTIATE_FREE_FUNCTIONS_REAL_ONLY(double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && (defined(ENABLE_COMPLEX_KERNELS) || defined(ENABLE_COMPLEX_BASIS_FUNCTIONS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>, std::complex<double>);
#endif

} // namespace Bempp
