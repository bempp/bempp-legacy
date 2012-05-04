// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_standard_local_assembler_factory_for_operators_on_surfaces_hpp
#define fiber_standard_local_assembler_factory_for_operators_on_surfaces_hpp

#include "local_assembler_factory.hpp"

#include "accuracy_options.hpp"
#include "opencl_options.hpp"
#include "standard_local_assembler_for_identity_operator_on_surface.hpp"
#include "standard_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "standard_local_assembler_for_grid_functions_on_surfaces.hpp"
#include "standard_evaluator_for_integral_operators.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename BasisFunctionType, typename ResultType,
typename GeometryFactory, typename Enable>
class StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase :
        public LocalAssemblerFactory<BasisFunctionType, ResultType,
        GeometryFactory, Enable>
{   
public:
    typedef LocalAssemblerFactory<BasisFunctionType, ResultType,
    GeometryFactory, Enable> Base;
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase() {
    }

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase(
            const AccuracyOptions& accuracyOptions) :
        m_accuracyOptions(accuracyOptions) {
    }

public:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIdentityOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        typedef StandardLocalAssemblerForIdentityOperatorOnSurface<
                BasisFunctionType, ResultType, GeometryFactory>
                LocalAssemblerForIdentityOperator_;
        return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                    new LocalAssemblerForIdentityOperator_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, trialExpression,
                        openClHandler));
    }

private:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperatorsImplRealKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Kernel<CoordinateType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const {
        typedef CoordinateType KernelType;
        typedef StandardLocalAssemblerForIntegralOperatorsOnSurfaces<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            LocalAssemblerForIntegralOperators_;
        return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                    new LocalAssemblerForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, kernel, trialExpression,
                        openClHandler, parallelisationOptions,
                        cacheSingularIntegrals,
                        this->accuracyOptions()));
    }

    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
    makeAssemblerForGridFunctionsImplRealUserFunction(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const Expression<CoordinateType>& testExpression,
            const Function<CoordinateType>& function,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        typedef CoordinateType UserFunctionType;
        typedef StandardLocalAssemblerForGridFunctionsOnSurfaces<
            BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>
            LocalAssemblerForGridFunctions_;
        return std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >(
                    new LocalAssemblerForGridFunctions_(
                        geometryFactory, rawGeometry,
                        testBases,
                        testExpression, function,
                        openClHandler));
    }

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperatorsImplRealKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<CoordinateType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        typedef CoordinateType KernelType;
        typedef StandardEvaluatorForIntegralOperators<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            EvaluatorForIntegralOperators_;
        return std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >(
                    new EvaluatorForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        trialBases,
                        kernel, trialExpression, argumentLocalCoefficients,
                        openClHandler,
                        this->accuracyOptions().singleRegular));
    }

public:
    const AccuracyOptions& accuracyOptions() const {
        return m_accuracyOptions;
    }

private:
    AccuracyOptions m_accuracyOptions;
};

// Complex ResultType
template <typename BasisFunctionType, typename ResultType,
typename GeometryFactory, typename Enable = void>
class StandardLocalAssemblerFactoryForOperatorsOnSurfaces :
        public StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase<
        BasisFunctionType, ResultType, GeometryFactory, Enable>
{
    typedef StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase<
    BasisFunctionType, ResultType, GeometryFactory, Enable> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces() : Base() {
    }

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit StandardLocalAssemblerFactoryForOperatorsOnSurfaces(
            const AccuracyOptions& accuracyOptions) :
        Base(accuracyOptions) {
    }

private:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperatorsImplComplexKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Kernel<ResultType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const {
        typedef ResultType KernelType;
        typedef StandardLocalAssemblerForIntegralOperatorsOnSurfaces<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            LocalAssemblerForIntegralOperators_;
        return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                    new LocalAssemblerForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, kernel, trialExpression,
                        openClHandler, parallelisationOptions,
                        cacheSingularIntegrals,
                        this->accuracyOptions()));
    }

    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
    makeAssemblerForGridFunctionsImplComplexUserFunction(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const Expression<CoordinateType>& testExpression,
            const Function<ResultType>& function,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        typedef ResultType UserFunctionType;
        typedef StandardLocalAssemblerForGridFunctionsOnSurfaces<
            BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>
            LocalAssemblerForGridFunctions_;
        return std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >(
                    new LocalAssemblerForGridFunctions_(
                        geometryFactory, rawGeometry,
                        testBases,
                        testExpression, function,
                        openClHandler));
    }

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperatorsImplComplexKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        typedef ResultType KernelType;
        typedef StandardEvaluatorForIntegralOperators<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
            EvaluatorForIntegralOperators_;
        return std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >(
                    new EvaluatorForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        trialBases,
                        kernel, trialExpression, argumentLocalCoefficients,
                        openClHandler,
                        this->accuracyOptions().singleRegular));
    }
};

// RealResultType
template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
class StandardLocalAssemblerFactoryForOperatorsOnSurfaces<
    BasisFunctionType, ResultType, GeometryFactory,
    typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type> :
        public StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase<
            BasisFunctionType, ResultType, GeometryFactory,
            typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type>
{
    typedef typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type Enable;
    typedef StandardLocalAssemblerFactoryForOperatorsOnSurfacesBase<
    BasisFunctionType, ResultType, GeometryFactory, Enable> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces() : Base() {
    }

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit StandardLocalAssemblerFactoryForOperatorsOnSurfaces(
            const AccuracyOptions& accuracyOptions) :
        Base(accuracyOptions) {
    }
};

} // namespace Fiber

#endif
