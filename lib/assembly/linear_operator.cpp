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
    m_testSpace(other.m_testSpace), m_trialSpace(other.m_trialSpace)
{
}

template <typename BasisFunctionType, typename ResultType>
LinearOperator<BasisFunctionType, ResultType>::~LinearOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
void LinearOperator<BasisFunctionType, ResultType>::assembleWeakForm(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry)
{
    m_weakForm = this->assembleDetachedWeakFormImpl(factory, options, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperator<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    return this->assembleDetachedWeakFormImpl(factory, options, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
bool LinearOperator<BasisFunctionType, ResultType>::isWeakFormAssembled() const
{
    return m_weakForm.get() != 0;
}

template <typename BasisFunctionType, typename ResultType>
const DiscreteLinearOperator<ResultType>&
LinearOperator<BasisFunctionType, ResultType>::weakForm() const
{
    if (!isWeakFormAssembled())
        throw std::runtime_error("LinearOperator::weakForm(): "
                                 "the weak form is not assembled");
    return *m_weakForm;
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperator<BasisFunctionType, ResultType>::detachWeakForm()
{
    return std::auto_ptr<DiscreteLinearOperator<ResultType> >(m_weakForm.release());
}

template <typename BasisFunctionType, typename ResultType>
void LinearOperator<BasisFunctionType, ResultType>::apply(
        const TranspositionMode trans,
        const GridFunction<BasisFunctionType, ResultType>& x_in,
        GridFunction<BasisFunctionType, ResultType>& y_inout,
        ResultType alpha, ResultType beta) const
{
    if (!isWeakFormAssembled())
        throw std::runtime_error("LinearOperator::apply(): "
                                 "the weak form is not assembled");

    // Sanity test
    if (&m_trialSpace != &x_in.space() || &m_testSpace != &y_inout.space())
        throw std::runtime_error("LinearOperator::apply(): Spaces don't match");

    // Extract coefficient vectors
    arma::Col<ResultType> xVals = x_in.coefficients();
    arma::Col<ResultType> yVals = y_inout.projections();

    // Apply operator and assign the result to y_inout's projections
    m_weakForm->apply(trans, xVals, yVals, alpha, beta);
    // TODO: make interfaces to the Trilinos and fallback
    // DiscreteLinearOperator::apply() compatible.
    // Perhaps by declaring an asPtrToBaseVector method in Vector...
    y_inout.setProjections(yVals);
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
LinearOperator<BasisFunctionType, ResultType>::testSpace() const
{
    return m_testSpace;
}

template <typename BasisFunctionType, typename ResultType>
const Space<BasisFunctionType>&
LinearOperator<BasisFunctionType, ResultType>::trialSpace() const
{
    return m_trialSpace;
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
    arma::Col<ResultType> coefficients(space.globalDofCount());
    coefficients.fill(0.);
    arma::Col<ResultType> projections(space.globalDofCount());
    projections.fill(0.);
    GridFunction<BasisFunctionType, ResultType> result(space, coefficients, projections);
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
