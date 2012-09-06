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

#ifndef fiber_numerical_quadrature_strategy_hpp
#define fiber_numerical_quadrature_strategy_hpp

#include "../common/common.hpp"

#include "quadrature_strategy.hpp"

#include "accuracy_options.hpp"

namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
typename GeometryFactory, typename Enable>
class NumericalQuadratureStrategyBase :
        public QuadratureStrategy<BasisFunctionType, ResultType,
        GeometryFactory, Enable>
{   
public:
    typedef QuadratureStrategy<BasisFunctionType, ResultType,
    GeometryFactory, Enable> Base;
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    NumericalQuadratureStrategyBase();

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit NumericalQuadratureStrategyBase(
            const AccuracyOptions& accuracyOptions);

public:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIdentityOperators(
            const shared_ptr<const GeometryFactory>& geometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
            const shared_ptr<const OpenClHandler>& openClHandler) const;

private:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperatorsImplRealKernel(
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
            const shared_ptr<const CollectionOfKernels<CoordinateType> >& kernels,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
            const shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, CoordinateType, ResultType> >& integral,
            const shared_ptr<const OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            bool cacheSingularIntegrals) const;

    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
    makeAssemblerForGridFunctionsImplRealUserFunction(
            const shared_ptr<const GeometryFactory>& geometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
            const shared_ptr<const Function<CoordinateType> >& function,
            const shared_ptr<const OpenClHandler>& openClHandler) const;

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperatorsImplRealKernel(
            const shared_ptr<const GeometryFactory>& geometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const CollectionOfKernels<CoordinateType> >& kernels,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
            const shared_ptr<const KernelTrialIntegral<BasisFunctionType, CoordinateType, ResultType> >& integral,
            const shared_ptr<const std::vector<std::vector<ResultType> > >& argumentLocalCoefficients,
            const shared_ptr<const OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions) const;

public:
    const AccuracyOptions& accuracyOptions() const;

private:
    AccuracyOptions m_accuracyOptions;
};

// Complex ResultType
template <typename BasisFunctionType, typename ResultType,
typename GeometryFactory, typename Enable = void>
class NumericalQuadratureStrategy :
        public NumericalQuadratureStrategyBase<
        BasisFunctionType, ResultType, GeometryFactory, Enable>
{
    typedef NumericalQuadratureStrategyBase<
    BasisFunctionType, ResultType, GeometryFactory, Enable> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    NumericalQuadratureStrategy();

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit NumericalQuadratureStrategy(
            const AccuracyOptions& accuracyOptions);

private:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperatorsImplComplexKernel(
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
            const shared_ptr<const CollectionOfKernels<ResultType> >& kernels,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
            const shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, ResultType, ResultType> >& integral,
            const shared_ptr<const OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            bool cacheSingularIntegrals) const;

    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
    makeAssemblerForGridFunctionsImplComplexUserFunction(
            const shared_ptr<const GeometryFactory>& geometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
            const shared_ptr<const Function<ResultType> >& function,
            const shared_ptr<const OpenClHandler>& openClHandler) const;

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperatorsImplComplexKernel(
            const shared_ptr<const GeometryFactory>& geometryFactory,
            const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
            const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const CollectionOfKernels<ResultType> >& kernels,
            const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
            const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType, ResultType> >& integral,
            const shared_ptr<const std::vector<std::vector<ResultType> > >& argumentLocalCoefficients,
            const shared_ptr<const OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions) const;
};

// RealResultType
template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
class NumericalQuadratureStrategy<
    BasisFunctionType, ResultType, GeometryFactory,
    typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type> :
        public NumericalQuadratureStrategyBase<
            BasisFunctionType, ResultType, GeometryFactory,
            typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type>
{
    typedef typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type Enable;
    typedef NumericalQuadratureStrategyBase<
    BasisFunctionType, ResultType, GeometryFactory, Enable> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    NumericalQuadratureStrategy();

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit NumericalQuadratureStrategy(
            const AccuracyOptions& accuracyOptions);
};

} // namespace Fiber

#endif
