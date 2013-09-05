// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_numerical_quadrature_strategy_imp_hpp
#define fiber_numerical_quadrature_strategy_imp_hpp

#include "numerical_quadrature_strategy.hpp"

#include "default_local_assembler_for_grid_functions_on_surfaces.hpp"
#include "default_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "default_local_assembler_for_local_operators_on_surfaces.hpp"
#include "default_local_assembler_for_potential_operators_on_surfaces.hpp"
#include "default_evaluator_for_integral_operators.hpp"
#include "default_quadrature_descriptor_selector_factory.hpp"

#include "default_double_quadrature_rule_family.hpp"
#include "default_single_quadrature_rule_family.hpp"

#include "default_test_trial_integral_imp.hpp"
#include "simple_test_trial_integrand_functor.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
NumericalQuadratureStrategyBase(
    const shared_ptr<const QuadratureDescriptorSelectorFactory<BasisFunctionType> >&
    quadratureDescriptorSelectorFactory,
    const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType> >&
    singleQuadratureRuleFamily,
    const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType> >&
    doubleQuadratureRuleFamily) :
    m_quadratureDescriptorSelectorFactory(
        quadratureDescriptorSelectorFactory),
    m_singleQuadratureRuleFamily(singleQuadratureRuleFamily),
    m_doubleQuadratureRuleFamily(doubleQuadratureRuleFamily)
{
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForOperators<ResultType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForIdentityOperators(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const OpenClHandler>& openClHandler) const
{
    typedef Fiber::SimpleTestTrialIntegrandFunctor<BasisFunctionType, ResultType>
            IntegrandFunctor;
    shared_ptr<TestTrialIntegral<BasisFunctionType, ResultType> > integral(
                new Fiber::DefaultTestTrialIntegral<IntegrandFunctor>(
                    IntegrandFunctor()));
    return makeAssemblerForLocalOperators(
                geometryFactory, rawGeometry,
                testShapesets, trialShapesets,
                testTransformations, trialTransformations, integral,
                openClHandler);
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForOperators<ResultType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForLocalOperators(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const TestTrialIntegral<BasisFunctionType, ResultType> >& integral,
        const shared_ptr<const OpenClHandler>& openClHandler) const
{
    typedef DefaultLocalAssemblerForLocalOperatorsOnSurfaces<
            BasisFunctionType, ResultType, GeometryFactory>
            LocalAssemblerForLocalOperators_;
    return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                new LocalAssemblerForLocalOperators_(
                    geometryFactory, rawGeometry,
                    testShapesets, trialShapesets,
                    testTransformations, trialTransformations, integral,
                    openClHandler,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForLocalOperators(
                        rawGeometry, testShapesets, trialShapesets),
                    this->singleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForOperators<ResultType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForIntegralOperatorsImplRealKernel(
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const CollectionOfKernels<CoordinateType> >& kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, CoordinateType, ResultType> >& integral,
        const shared_ptr<const OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        VerbosityLevel::Level verbosityLevel,
        bool cacheSingularIntegrals) const
{
    typedef CoordinateType KernelType;
    typedef DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            LocalAssemblerForIntegralOperators_;
    return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                new LocalAssemblerForIntegralOperators_(
                    testGeometryFactory, trialGeometryFactory,
                    testRawGeometry, trialRawGeometry,
                    testShapesets, trialShapesets,
                    testTransformations, kernels, trialTransformations, integral,
                    openClHandler, parallelizationOptions,
                    verbosityLevel,
                    cacheSingularIntegrals,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForIntegralOperators(
                        testRawGeometry, trialRawGeometry,
                        testShapesets, trialShapesets),
                    this->doubleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForGridFunctionsImplRealUserFunction(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const Function<CoordinateType> >& function,
        const shared_ptr<const OpenClHandler>& openClHandler) const
{
    typedef CoordinateType UserFunctionType;
    typedef DefaultLocalAssemblerForGridFunctionsOnSurfaces<
            BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>
            LocalAssemblerForGridFunctions_;
    return std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >(
                new LocalAssemblerForGridFunctions_(
                    geometryFactory, rawGeometry,
                    testShapesets,
                    testTransformations, function,
                    openClHandler,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForGridFunctions(
                        rawGeometry, testShapesets),
                    this->singleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeEvaluatorForIntegralOperatorsImplRealKernel(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfKernels<CoordinateType> >& kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const KernelTrialIntegral<BasisFunctionType, CoordinateType, ResultType> >& integral,
        const shared_ptr<const std::vector<std::vector<ResultType> > >& argumentLocalCoefficients,
        const shared_ptr<const OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions) const
{
    typedef CoordinateType KernelType;
    typedef DefaultEvaluatorForIntegralOperators<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            EvaluatorForIntegralOperators_;
    return std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >(
                new EvaluatorForIntegralOperators_(
                    geometryFactory, rawGeometry,
                    trialShapesets,
                    kernels, trialTransformations, integral,
                    argumentLocalCoefficients,
                    openClHandler,
                    parallelizationOptions,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForPotentialOperators(
                        rawGeometry, trialShapesets),
                    this->singleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForPotentialOperators<ResultType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForPotentialOperatorsImplRealKernel(
        const arma::Mat<CoordinateType>& evaluationPoints,
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfKernels<CoordinateType> >& kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const KernelTrialIntegral<BasisFunctionType, CoordinateType, ResultType> >& integral,
        const shared_ptr<const OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        VerbosityLevel::Level verbosityLevel) const
{
    typedef CoordinateType KernelType;
    typedef DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            LocalAssemblerForPotentialOperators_;
    return std::auto_ptr<LocalAssemblerForPotentialOperators<ResultType> >(
                new LocalAssemblerForPotentialOperators_(
                    evaluationPoints,
                    geometryFactory, rawGeometry,
                    trialShapesets,
                    kernels, trialTransformations, integral,
                    parallelizationOptions,
                    verbosityLevel,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForPotentialOperators(
                        rawGeometry, trialShapesets),
                    this->singleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
shared_ptr<const QuadratureDescriptorSelectorFactory<BasisFunctionType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
quadratureDescriptorSelectorFactory() const
{
    return m_quadratureDescriptorSelectorFactory;
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
shared_ptr<const DoubleQuadratureRuleFamily<
typename NumericalQuadratureStrategyBase<
    BasisFunctionType, ResultType, GeometryFactory, Enable>::CoordinateType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
doubleQuadratureRuleFamily() const
{
    return m_doubleQuadratureRuleFamily;
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
shared_ptr<const SingleQuadratureRuleFamily<
typename NumericalQuadratureStrategyBase<
    BasisFunctionType, ResultType, GeometryFactory, Enable>::CoordinateType> >
NumericalQuadratureStrategyBase<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
singleQuadratureRuleFamily() const
{
    return m_singleQuadratureRuleFamily;
}

// Complex ResultType
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
NumericalQuadratureStrategy() :
    Base(
        boost::make_shared<
            DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType> >(),
        boost::make_shared<
            DefaultSingleQuadratureRuleFamily<CoordinateType> >(),
        boost::make_shared<
            DefaultDoubleQuadratureRuleFamily<CoordinateType> >())
{
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
NumericalQuadratureStrategy(
        const AccuracyOptionsEx& accuracyOptions) :
    Base(
        boost::make_shared<
            DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType> >(
                accuracyOptions),
        boost::make_shared<
            DefaultSingleQuadratureRuleFamily<CoordinateType> >(),
        boost::make_shared<
            DefaultDoubleQuadratureRuleFamily<CoordinateType> >())
{
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
NumericalQuadratureStrategy(
        const AccuracyOptionsEx& accuracyOptions,
        const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType> >&
        singleQuadratureRuleFamily,
        const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType> >&
        doubleQuadratureRuleFamily) :
    Base(
        boost::make_shared<
            DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType> >(accuracyOptions),
        singleQuadratureRuleFamily,
        doubleQuadratureRuleFamily)
{
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
NumericalQuadratureStrategy(
    const shared_ptr<const QuadratureDescriptorSelectorFactory<BasisFunctionType> >&
    quadratureDescriptorSelectorFactory,
    const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType> >&
    singleQuadratureRuleFamily,
    const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType> >&
    doubleQuadratureRuleFamily) :
    Base(
        quadratureDescriptorSelectorFactory,
        singleQuadratureRuleFamily,
        doubleQuadratureRuleFamily)
{
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForOperators<ResultType> >
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForIntegralOperatorsImplComplexKernel(
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const CollectionOfKernels<ResultType> >& kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, ResultType, ResultType> >& integral,
        const shared_ptr<const OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        VerbosityLevel::Level verbosityLevel,
        bool cacheSingularIntegrals) const
{
    typedef ResultType KernelType;
    typedef DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            LocalAssemblerForIntegralOperators_;
    return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                new LocalAssemblerForIntegralOperators_(
                    testGeometryFactory, trialGeometryFactory,
                    testRawGeometry, trialRawGeometry,
                    testShapesets, trialShapesets,
                    testTransformations, kernels, trialTransformations, integral,
                    openClHandler, parallelizationOptions,
                    verbosityLevel,
                    cacheSingularIntegrals,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForIntegralOperators(
                        testRawGeometry, trialRawGeometry,
                        testShapesets, trialShapesets),
                    this->doubleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForGridFunctionsImplComplexUserFunction(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const Function<ResultType> >& function,
        const shared_ptr<const OpenClHandler>& openClHandler) const
{
    typedef ResultType UserFunctionType;
    typedef DefaultLocalAssemblerForGridFunctionsOnSurfaces<
            BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>
            LocalAssemblerForGridFunctions_;
    return std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >(
                new LocalAssemblerForGridFunctions_(
                    geometryFactory, rawGeometry,
                    testShapesets,
                    testTransformations, function,
                    openClHandler,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForGridFunctions(
                        rawGeometry, testShapesets),
                    this->singleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeEvaluatorForIntegralOperatorsImplComplexKernel(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfKernels<ResultType> >& kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType, ResultType> >& integral,
        const shared_ptr<const std::vector<std::vector<ResultType> > >& argumentLocalCoefficients,
        const shared_ptr<const OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions) const
{
    typedef ResultType KernelType;
    typedef DefaultEvaluatorForIntegralOperators<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            EvaluatorForIntegralOperators_;
    return std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >(
                new EvaluatorForIntegralOperators_(
                    geometryFactory, rawGeometry,
                    trialShapesets,
                    kernels, trialTransformations, integral,
                    argumentLocalCoefficients,
                    openClHandler,
                    parallelizationOptions,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForPotentialOperators(
                        rawGeometry, trialShapesets),
                    this->singleQuadratureRuleFamily()));
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
std::auto_ptr<LocalAssemblerForPotentialOperators<ResultType> >
NumericalQuadratureStrategy<
BasisFunctionType, ResultType, GeometryFactory, Enable>::
makeAssemblerForPotentialOperatorsImplComplexKernel(
        const arma::Mat<CoordinateType>& evaluationPoints,
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const CollectionOfKernels<ResultType> >& kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType, ResultType> >& integral,
        const shared_ptr<const OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        VerbosityLevel::Level verbosityLevel) const
{
    typedef ResultType KernelType;
    typedef DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            LocalAssemblerForPotentialOperators_;
    return std::auto_ptr<LocalAssemblerForPotentialOperators<ResultType> >(
                new LocalAssemblerForPotentialOperators_(
                    evaluationPoints,
                    geometryFactory, rawGeometry,
                    trialShapesets,
                    kernels, trialTransformations, integral,
                    parallelizationOptions,
                    verbosityLevel,
                    this->quadratureDescriptorSelectorFactory()->
                    makeQuadratureDescriptorSelectorForPotentialOperators(
                        rawGeometry, trialShapesets),
                    this->singleQuadratureRuleFamily()));
}

// template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
// NumericalQuadratureStrategy<
//         BasisFunctionType, ResultType, GeometryFactory,
//         typename boost::enable_if<
//         boost::is_same<ResultType,
//         typename ScalarTraits<ResultType>::RealType> >::type>::
// NumericalQuadratureStrategy() :
//     Base()
// {
// }

// template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
// NumericalQuadratureStrategy<
//         BasisFunctionType, ResultType, GeometryFactory,
//         typename boost::enable_if<
//         boost::is_same<ResultType,
//         typename ScalarTraits<ResultType>::RealType> >::type>::
// NumericalQuadratureStrategy(
//         const AccuracyOptionsEx& accuracyOptions) :
//     Base(accuracyOptions)
// {
// }

template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
NumericalQuadratureStrategy<
        BasisFunctionType, ResultType, GeometryFactory,
        typename boost::enable_if<
        boost::is_same<ResultType,
        typename ScalarTraits<ResultType>::RealType> >::type>::
NumericalQuadratureStrategy() :
    Base(
        boost::make_shared<
            DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType> >(),
        boost::make_shared<
            DefaultSingleQuadratureRuleFamily<CoordinateType> >(),
        boost::make_shared<
            DefaultDoubleQuadratureRuleFamily<CoordinateType> >())
{
}

template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
NumericalQuadratureStrategy<
        BasisFunctionType, ResultType, GeometryFactory,
        typename boost::enable_if<
        boost::is_same<ResultType,
        typename ScalarTraits<ResultType>::RealType> >::type>::
NumericalQuadratureStrategy(
        const AccuracyOptionsEx& accuracyOptions) :
    Base(
        boost::make_shared<
            DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType> >(
                accuracyOptions),
        boost::make_shared<
            DefaultSingleQuadratureRuleFamily<CoordinateType> >(),
        boost::make_shared<
            DefaultDoubleQuadratureRuleFamily<CoordinateType> >())
{
}

template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
NumericalQuadratureStrategy<
        BasisFunctionType, ResultType, GeometryFactory,
        typename boost::enable_if<
        boost::is_same<ResultType,
        typename ScalarTraits<ResultType>::RealType> >::type>::
NumericalQuadratureStrategy(
        const AccuracyOptionsEx& accuracyOptions,
        const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType> >&
        singleQuadratureRuleFamily,
        const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType> >&
        doubleQuadratureRuleFamily) :
    Base(
        boost::make_shared<
            DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType> >(accuracyOptions),
        singleQuadratureRuleFamily,
        doubleQuadratureRuleFamily)
{
}

template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
NumericalQuadratureStrategy<
        BasisFunctionType, ResultType, GeometryFactory,
        typename boost::enable_if<
        boost::is_same<ResultType,
        typename ScalarTraits<ResultType>::RealType> >::type>::
NumericalQuadratureStrategy(
    const shared_ptr<const QuadratureDescriptorSelectorFactory<BasisFunctionType> >&
    quadratureDescriptorSelectorFactory,
    const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType> >&
    singleQuadratureRuleFamily,
    const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType> >&
    doubleQuadratureRuleFamily) :
    Base(
        quadratureDescriptorSelectorFactory,
        singleQuadratureRuleFamily,
        doubleQuadratureRuleFamily)
{
}

} // namespace Fiber

#endif
