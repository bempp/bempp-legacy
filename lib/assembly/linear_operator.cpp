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

#include "linear_operator.hpp"
#include "linear_operator_superposition.hpp"
#include "discrete_linear_operator.hpp"
#include "grid_function.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include <stdexcept>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
LinearOperator<BasisFunctionType, ResultType>::
LinearOperator(const Space<BasisFunctionType>& testSpace,
               const Space<BasisFunctionType>& trialSpace) :
    m_testSpace(testSpace), m_trialSpace(trialSpace)
{
}

template <typename BasisFunctionType, typename ResultType>
LinearOperator<BasisFunctionType, ResultType>::LinearOperator(
        const LinearOperator<BasisFunctionType, ResultType>& other) :
    m_testSpace(other.m_testSpace), m_trialSpace(other.m_trialSpace),
    m_localOperators(other.m_localOperators), m_multipliers(other.m_multipliers)
{
}

template <typename BasisFunctionType, typename ResultType>
LinearOperator<BasisFunctionType, ResultType>::~LinearOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
void LinearOperator<BasisFunctionType, ResultType>::assemble(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options)
{
    m_discreteOperator = this->assembleWeakForm(factory, options);
}

template <typename BasisFunctionType, typename ResultType>
bool LinearOperator<BasisFunctionType, ResultType>::isAssembled() const
{
    return m_discreteOperator.get() != NULL;
}

template <typename BasisFunctionType, typename ResultType>
void LinearOperator<BasisFunctionType, ResultType>::apply(
        const TranspositionMode trans,
        const GridFunction<BasisFunctionType, ResultType>& x_in,
        GridFunction<BasisFunctionType, ResultType>& y_inout,
        ResultType alpha, ResultType beta) const
{
    if (!this->isAssembled())
        throw std::runtime_error("LinearOperator::apply(): "
                                 "operator is not assembled");

    // Sanity test
    if (&m_trialSpace != &x_in.space() || &m_testSpace != &y_inout.space())
        throw std::runtime_error("LinearOperator::apply(): Spaces don't match");

    // Extract coefficient vectors
    arma::Col<ResultType> xVals = x_in.coefficients().asArmadilloVector();
    arma::Col<ResultType> yVals = y_inout.coefficients().asArmadilloVector();

    // Apply operator and assign the result to y_inout's coefficients
    m_discreteOperator->apply(trans, xVals, yVals, alpha, beta);
    // Note: all these armadillo<->vector conversions are horribly
    // inefficient and unnecessary. But in order to get rid of them, we need
    // first to make interfaces to the Trilinos and fallback
    // DiscreteLinearOperator::apply() compatible.
    // Perhaps by declaring an asPtrToBaseVector method in Vector...
    y_inout.setCoefficients(Vector<ResultType>(yVals));
}

template <typename BasisFunctionType, typename ResultType>
const DiscreteLinearOperator<ResultType>&
LinearOperator<BasisFunctionType, ResultType>::assembledDiscreteLinearOperator() const
{
    if (!isAssembled())
        throw std::runtime_error("LinearOperator::assembledDiscreteLinearOperator(): "
                                 "operator is not assembled");
    return *m_discreteOperator;
}

template <typename BasisFunctionType, typename ResultType>
const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>&
LinearOperator<BasisFunctionType, ResultType>::localOperators() const
{
    return m_localOperators;
}

template <typename BasisFunctionType, typename ResultType>
const std::vector<ResultType>& LinearOperator<BasisFunctionType, ResultType>::multipliers() const
{
    return m_multipliers;
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>& LinearOperator<BasisFunctionType, ResultType>::testSpace() const
{
    return m_testSpace;
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>& LinearOperator<BasisFunctionType, ResultType>::trialSpace() const
{
    return m_trialSpace;
//    LinearOperatorSuperposition<BasisFunctionType, ResultType> ls =
//            *this * static_cast<ResultType>(5.);
}

template <typename BasisFunctionType, typename ResultType>
void LinearOperator<BasisFunctionType, ResultType>::addLocalOperatorsAndMultipliers(
        const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>&
        localOperators,
        const std::vector<ResultType>& multipliers)
{
    m_localOperators.insert(m_localOperators.end(),
                            localOperators.begin(), localOperators.end());
    m_multipliers.insert(m_multipliers.end(),
                         multipliers.begin(), multipliers.end());
}

template <typename BasisFunctionType, typename ResultType>
void
LinearOperator<BasisFunctionType, ResultType>::collectDataForAssemblerConstruction(
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
    Helper::collectGridData(m_testSpace.grid(),
                            testRawGeometry, testGeometryFactory);
    if (&m_testSpace.grid() == &m_trialSpace.grid()) {
        trialRawGeometry = testRawGeometry;
        trialGeometryFactory = testGeometryFactory;
    } else
        Helper::collectGridData(m_trialSpace.grid(),
                                trialRawGeometry, trialGeometryFactory);

    // Construct the OpenClHandler
    Helper::makeOpenClHandler(options.parallelisationOptions().openClOptions(),
                              testRawGeometry, trialRawGeometry, openClHandler);

    // Get pointers to test and trial bases of each element
    Helper::collectBases(m_testSpace, testBases);
    if (&m_testSpace == &m_trialSpace)
        trialBases = testBases;
    else
        Helper::collectBases(m_trialSpace, trialBases);

    cacheSingularIntegrals =
            (options.singularIntegralCaching() == AssemblyOptions::YES ||
             (options.singularIntegralCaching() == AssemblyOptions::AUTO));
}

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator+(
        const LinearOperator<BasisFunctionType, ResultType>& op1,
        const LinearOperator<BasisFunctionType, ResultType>& op2)
{
    return LinearOperatorSuperposition<BasisFunctionType, ResultType>(op1, op2);
}

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator-(
        const LinearOperator<BasisFunctionType, ResultType>& op1,
        const LinearOperator<BasisFunctionType, ResultType>& op2)
{
    return op1 + (static_cast<ResultType>(-1.) * op2);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>, std::complex<double> >,
        ScalarType>,
    LinearOperatorSuperposition<BasisFunctionType, ResultType> >::type
operator*(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar)
{
    return LinearOperatorSuperposition<BasisFunctionType, ResultType>(
                op, static_cast<ResultType>(scalar));
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator*(
        const ScalarType& scalar,
        const LinearOperator<BasisFunctionType, ResultType>& op)
{
     return operator*(op, scalar);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
LinearOperatorSuperposition<BasisFunctionType, ResultType> operator/(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const ScalarType& scalar)
{
    if (scalar == static_cast<ScalarType>(0.))
        throw std::runtime_error("LinearOperatorSuperposition::operator/(): "
                                 "Division by zero");

    return LinearOperatorSuperposition<BasisFunctionType, ResultType>(
                op, static_cast<ResultType>(static_cast<ScalarType>(1.) / scalar));
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> operator*(
        const LinearOperator<BasisFunctionType, ResultType>& op,
        const GridFunction<BasisFunctionType, ResultType>& fun)
{
    const Space<BasisFunctionType>& space = op.testSpace();
    arma::Col<ResultType> coeffs(space.globalDofCount());
    coeffs.fill(0.);
    GridFunction<BasisFunctionType, ResultType> result(space, coeffs);
    op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
    return result;
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT) \
    template LinearOperatorSuperposition<BASIS, RESULT> operator+( \
    const LinearOperator<BASIS, RESULT>& op1, \
    const LinearOperator<BASIS, RESULT>& op2); \
    template LinearOperatorSuperposition<BASIS, RESULT> operator-( \
    const LinearOperator<BASIS, RESULT>& op1, \
    const LinearOperator<BASIS, RESULT>& op2); \
    template GridFunction<BASIS, RESULT> operator*( \
    const LinearOperator<BASIS, RESULT>& op, \
    const GridFunction<BASIS, RESULT>& fun)
#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(BASIS, RESULT, SCALAR) \
    template LinearOperatorSuperposition<BASIS, RESULT> operator*( \
    const LinearOperator<BASIS, RESULT>& op, const SCALAR& scalar); \
    template LinearOperatorSuperposition<BASIS, RESULT> operator*( \
    const SCALAR& scalar, const LinearOperator<BASIS, RESULT>& op); \
    template LinearOperatorSuperposition<BASIS, RESULT> operator/( \
    const LinearOperator<BASIS, RESULT>& op, const SCALAR& scalar)

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);

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
