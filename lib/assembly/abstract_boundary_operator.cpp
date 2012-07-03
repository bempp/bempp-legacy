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

#include "abstract_boundary_operator.hpp"
#include "discrete_boundary_operator.hpp"
#include "grid_function.hpp"
// #include "boundary_operator_composition.hpp"
#include "abstract_boundary_operator_sum.hpp"
#include "local_assembler_construction_helper.hpp"
#include "scaled_abstract_boundary_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include <stdexcept>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperator<BasisFunctionType, ResultType>::
AbstractBoundaryOperator(const Space<BasisFunctionType>& domain,
               const Space<BasisFunctionType>& range,
               const Space<BasisFunctionType>& dualToRange,
               const std::string& label) :
    m_domain(domain), m_range(range), m_dualToRange(dualToRange),
    m_label(label)
{
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperator<BasisFunctionType, ResultType>::~AbstractBoundaryOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
void AbstractBoundaryOperator<BasisFunctionType, ResultType>::assembleWeakForm(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry)
{
    m_weakForm = this->assembleWeakFormImpl(factory, options, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
bool AbstractBoundaryOperator<BasisFunctionType, ResultType>::isWeakFormAssembled() const
{
    return m_weakForm.get() != 0;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType> >
AbstractBoundaryOperator<BasisFunctionType, ResultType>::weakForm() const
{
    if (!isWeakFormAssembled())
        throw std::runtime_error("AbstractBoundaryOperator::weakForm(): "
                                 "the weak form is not assembled");
    return m_weakForm;
}

template <typename BasisFunctionType, typename ResultType>
void
AbstractBoundaryOperator<BasisFunctionType, ResultType>::resetWeakForm()
{
    m_weakForm.reset();
}

template <typename BasisFunctionType, typename ResultType>
void AbstractBoundaryOperator<BasisFunctionType, ResultType>::apply(
        const TranspositionMode trans,
        const GridFunction<BasisFunctionType, ResultType>& x_in,
        GridFunction<BasisFunctionType, ResultType>& y_inout,
        ResultType alpha, ResultType beta) const
{
    if (!isWeakFormAssembled())
        throw std::runtime_error("AbstractBoundaryOperator::apply(): "
                                 "the weak form is not assembled");

    // Sanity test
    if (&m_domain != &x_in.space() ||
            &m_range != &y_inout.space() ||
            &m_dualToRange != &y_inout.dualSpace())
        throw std::runtime_error("AbstractBoundaryOperator::apply(): Spaces don't match");

    // Extract coefficient vectors
    arma::Col<ResultType> xVals = x_in.coefficients();
    arma::Col<ResultType> yVals = y_inout.projections();

    // Apply operator and assign the result to y_inout's projections
    m_weakForm->apply(trans, xVals, yVals, alpha, beta);
    // TODO: make interfaces to the Trilinos and fallback
    // DiscreteBoundaryOperator::apply() compatible.
    // Perhaps by declaring an asPtrToBaseVector method in Vector...
    y_inout.setProjections(yVals);
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
AbstractBoundaryOperator<BasisFunctionType, ResultType>::domain() const
{
    return m_domain;
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
AbstractBoundaryOperator<BasisFunctionType, ResultType>::range() const
{
    return m_range;
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
AbstractBoundaryOperator<BasisFunctionType, ResultType>::dualToRange() const
{
    return m_dualToRange;
}

template <typename BasisFunctionType, typename ResultType>
std::string
AbstractBoundaryOperator<BasisFunctionType, ResultType>::label() const
{
    return m_label;
}

template <typename BasisFunctionType, typename ResultType>
void
AbstractBoundaryOperator<BasisFunctionType, ResultType>::setLabel(
        const std::string &newLabel)
{
    m_label = newLabel;
}

template <typename BasisFunctionType, typename ResultType>
void
AbstractBoundaryOperator<BasisFunctionType, ResultType>::collectDataForAssemblerConstruction(
        const AssemblyOptions& options,
        shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        shared_ptr<Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        shared_ptr<GeometryFactory>& testGeometryFactory,
        shared_ptr<GeometryFactory>& trialGeometryFactory,
        shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        shared_ptr<std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        shared_ptr<Fiber::OpenClHandler>& openClHandler,
        bool& cacheSingularIntegrals) const
{
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;
    typedef LocalAssemblerConstructionHelper Helper;

    // Collect grid data
    Helper::collectGridData(m_dualToRange.grid(),
                            testRawGeometry, testGeometryFactory);
    if (&m_dualToRange.grid() == &m_domain.grid()) {
        trialRawGeometry = testRawGeometry;
        trialGeometryFactory = testGeometryFactory;
    } else
        Helper::collectGridData(m_domain.grid(),
                                trialRawGeometry, trialGeometryFactory);

    // Construct the OpenClHandler
    Helper::makeOpenClHandler(options.parallelisationOptions().openClOptions(),
                              testRawGeometry, trialRawGeometry, openClHandler);

    // Get pointers to test and trial bases of each element
    Helper::collectBases(m_dualToRange, testBases);
    if (&m_dualToRange == &m_domain)
        trialBases = testBases;
    else
        Helper::collectBases(m_domain, trialBases);

    cacheSingularIntegrals =
            (options.singularIntegralCaching() == AssemblyOptions::YES ||
             (options.singularIntegralCaching() == AssemblyOptions::AUTO));
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> operator+(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op1,
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op2)
{
    return AbstractBoundaryOperatorSum<BasisFunctionType, ResultType>(op1, op2);
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> operator-(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op1,
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op2)
{
    return op1 + (static_cast<ResultType>(-1.) * op2);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> >::type
operator*(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar)
{
    return ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType>(
                static_cast<ResultType>(scalar), op);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar,
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op)
{
     return operator*(op, scalar);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> operator/(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar)
{
    if (scalar == static_cast<ScalarType>(0.))
        throw std::runtime_error("ScaledAbstractBoundaryOperator::operator/(): "
                                 "Division by zero");

    return ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType>(
                static_cast<ResultType>(static_cast<ScalarType>(1.) / scalar), op);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& op,
        const GridFunction<BasisFunctionType, ResultType>& fun)
{
    typedef GridFunction<BasisFunctionType, ResultType> GF;

    const Space<BasisFunctionType>& space = op.range();
    const Space<BasisFunctionType>& dualSpace = op.dualToRange();
    arma::Col<ResultType> coefficients(space.globalDofCount());
    coefficients.fill(0.);
    arma::Col<ResultType> projections(dualSpace.globalDofCount());
    projections.fill(0.);
//    GF result(space, dualSpace, coefficients, projections);
    GF result(space, dualSpace, projections, GF::PROJECTIONS);
//    GF result(space, dualSpace, coefficients, GF::COEFFICIENTS);
    op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
    return result;
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
    template AbstractBoundaryOperatorSum<BASIS, RESULT> operator+( \
    const AbstractBoundaryOperator<BASIS, RESULT>& op1, \
    const AbstractBoundaryOperator<BASIS, RESULT>& op2); \
    template AbstractBoundaryOperatorSum<BASIS, RESULT> operator-( \
    const AbstractBoundaryOperator<BASIS, RESULT>& op1, \
    const AbstractBoundaryOperator<BASIS, RESULT>& op2); \
    template GridFunction<BASIS, RESULT> operator*( \
    const AbstractBoundaryOperator<BASIS, RESULT>& op, \
    const GridFunction<BASIS, RESULT>& fun)
#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(BASIS, RESULT, SCALAR) \
    template ScaledAbstractBoundaryOperator<BASIS, RESULT> operator*( \
    const AbstractBoundaryOperator<BASIS, RESULT>& op, const SCALAR& scalar); \
    template ScaledAbstractBoundaryOperator<BASIS, RESULT> operator*( \
    const SCALAR& scalar, const AbstractBoundaryOperator<BASIS, RESULT>& op); \
    template ScaledAbstractBoundaryOperator<BASIS, RESULT> operator/( \
    const AbstractBoundaryOperator<BASIS, RESULT>& op, const SCALAR& scalar)

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

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
