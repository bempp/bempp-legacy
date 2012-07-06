#include "boundary_operator.hpp"

#include "abstract_boundary_operator_sum.hpp"
#include "discrete_boundary_operator.hpp"
#include "context.hpp"
#include "grid_function.hpp"
#include "scaled_abstract_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>::BoundaryOperator(
        const shared_ptr<const Context<
        BasisFunctionType, ResultType> >& context,
        const shared_ptr<const AbstractBoundaryOperator<
        BasisFunctionType, ResultType> >& abstractOp) :
    m_context(context), m_abstractOp(abstractOp)
{
    if (!m_context)
        throw std::invalid_argument("BoundaryOperator::BoundaryOperator(): "
                                    "context must not be null");
    if (!m_abstractOp)
        throw std::invalid_argument("BoundaryOperator::BoundaryOperator(): "
                                    "abstractOp must not be null");
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType> >
BoundaryOperator<BasisFunctionType, ResultType>::abstractOperator() const
{
    return m_abstractOp;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Context<BasisFunctionType, ResultType> >
BoundaryOperator<BasisFunctionType, ResultType>::context() const
{
    return m_context;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType> >
BoundaryOperator<BasisFunctionType, ResultType>::weakForm() const
{
    if (!m_weakForm.get()) {
        tbb::mutex::scoped_lock lock(m_weakFormMutex);
        if (!m_weakForm.get()) {
            m_weakForm = m_context->getWeakForm(*m_abstractOp);
            assert(m_weakForm);
        }
    }
    return m_weakForm;
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
BoundaryOperator<BasisFunctionType, ResultType>::domain() const
{
    return m_abstractOp->domain();
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
BoundaryOperator<BasisFunctionType, ResultType>::range() const
{
    return m_abstractOp->range();
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
BoundaryOperator<BasisFunctionType, ResultType>::dualToRange() const
{
    return m_abstractOp->dualToRange();
}

template <typename BasisFunctionType, typename ResultType>
std::string
BoundaryOperator<BasisFunctionType, ResultType>::label() const
{
    return m_abstractOp->label();
}

template <typename BasisFunctionType, typename ResultType>
void BoundaryOperator<BasisFunctionType, ResultType>::apply(
        const TranspositionMode trans,
        const GridFunction<BasisFunctionType, ResultType>& x_in,
        GridFunction<BasisFunctionType, ResultType>& y_inout,
        ResultType alpha, ResultType beta) const
{
    // Sanity test
    if (&m_abstractOp->domain() != &x_in.space() ||
            &m_abstractOp->range() != &y_inout.space() ||
            &m_abstractOp->dualToRange() != &y_inout.dualSpace())
        throw std::runtime_error("AbstractBoundaryOperator::apply(): Spaces don't match");

    // Extract coefficient vectors
    arma::Col<ResultType> xVals = x_in.coefficients();
    arma::Col<ResultType> yVals = y_inout.projections();

    // Apply operator and assign the result to y_inout's projections
    weakForm()->apply(trans, xVals, yVals, alpha, beta);
    // TODO: make interfaces to the Trilinos and fallback
    // DiscreteBoundaryOperator::apply() compatible.
    // Perhaps by declaring an asPtrToBaseVector method in Vector...
    y_inout.setProjections(yVals);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> operator+(
        const BoundaryOperator<BasisFunctionType, ResultType>& op1,
        const BoundaryOperator<BasisFunctionType, ResultType>& op2)
{
    typedef AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> Sum;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                op1.context(),
                boost::make_shared<Sum>(op1, op2));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> operator-(
        const BoundaryOperator<BasisFunctionType, ResultType>& op1,
        const BoundaryOperator<BasisFunctionType, ResultType>& op2)
{
    return op1 + (static_cast<ResultType>(-1.) * op2);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    BoundaryOperator<BasisFunctionType, ResultType> >::type
operator*(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar)
{
    typedef ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> ScaledOp;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                op.context(),
                boost::make_shared<ScaledOp>(static_cast<ResultType>(scalar), op));
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar,
        const BoundaryOperator<BasisFunctionType, ResultType>& op)
{
     return operator*(op, scalar);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType> operator/(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar)
{
    if (scalar == static_cast<ScalarType>(0.))
        throw std::runtime_error("ScaledAbstractBoundaryOperator::operator/(): "
                                 "Division by zero");
    return operator*(op, static_cast<ResultType>(static_cast<ScalarType>(1.) / scalar));
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        const GridFunction<BasisFunctionType, ResultType>& fun)
{
    typedef GridFunction<BasisFunctionType, ResultType> GF;

    const Space<BasisFunctionType>& space = op.abstractOperator()->range();
    const Space<BasisFunctionType>& dualSpace = op.abstractOperator()->dualToRange();
    arma::Col<ResultType> coefficients(space.globalDofCount());
    coefficients.fill(0.);
    arma::Col<ResultType> projections(dualSpace.globalDofCount());
    projections.fill(0.);
    GF result(op.context(), space, dualSpace, projections, GF::PROJECTIONS);
    op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
    return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> operator+( \
    const BoundaryOperator<BASIS, RESULT>& op1, \
    const BoundaryOperator<BASIS, RESULT>& op2); \
    template BoundaryOperator<BASIS, RESULT> operator-( \
    const BoundaryOperator<BASIS, RESULT>& op1, \
    const BoundaryOperator<BASIS, RESULT>& op2); \
    template GridFunction<BASIS, RESULT> operator*( \
    const BoundaryOperator<BASIS, RESULT>& op, \
    const GridFunction<BASIS, RESULT>& fun)
#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(BASIS, RESULT, SCALAR) \
    template BoundaryOperator<BASIS, RESULT> operator*( \
    const BoundaryOperator<BASIS, RESULT>& op, const SCALAR& scalar); \
    template BoundaryOperator<BASIS, RESULT> operator*( \
    const SCALAR& scalar, const BoundaryOperator<BASIS, RESULT>& op); \
    template BoundaryOperator<BASIS, RESULT> operator/( \
    const BoundaryOperator<BASIS, RESULT>& op, const SCALAR& scalar)

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(
        float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, float, double);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
INSTANTIATE_FREE_FUNCTIONS(
        float, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        float, std::complex<float>, std::complex<double>);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
INSTANTIATE_FREE_FUNCTIONS(
        std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<float>, std::complex<float>, std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(
        double, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, double, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, double, double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
INSTANTIATE_FREE_FUNCTIONS(
        double, std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        double, std::complex<double>, std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
INSTANTIATE_FREE_FUNCTIONS(
        std::complex<double>, std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(
        std::complex<double>, std::complex<double>, std::complex<double>);
#endif


} // namespace Bempp
