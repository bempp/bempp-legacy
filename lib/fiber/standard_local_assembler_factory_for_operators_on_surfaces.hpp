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

#include "accuracy_options.hpp"
#include "local_assembler_factory.hpp"
#include "opencl_options.hpp"
#include "standard_local_assembler_for_identity_operator_on_surface.hpp"
#include "standard_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "standard_local_assembler_for_grid_functions_on_surfaces.hpp"
#include "standard_evaluator_for_integral_operators.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename BasisValueType, typename ResultType,
typename GeometryFactory>
class StandardLocalAssemblerFactoryForOperatorsOnSurfaces :
        public LocalAssemblerFactory<BasisValueType, ResultType,
        GeometryFactory>
{   
public:
    typedef LocalAssemblerFactory<BasisValueType, ResultType,
    GeometryFactory> Base;
    typedef typename Base::CoordinateType CoordinateType;

    /** \brief Construct a local assembler factory with default accuracy settings. */
    StandardLocalAssemblerFactoryForOperatorsOnSurfaces() {
    }

    /** \brief Construct a local assembler factory with specified accuracy settings. */
    explicit StandardLocalAssemblerFactoryForOperatorsOnSurfaces(
            const AccuracyOptions& accuracyOptions) :
        m_accuracyOptions(accuracyOptions) {
    }

private:
    typedef StandardLocalAssemblerForIntegralOperatorsOnSurfaces<
        BasisValueType, ResultType, GeometryFactory>
        LocalAssemblerForIntegralOperators_;
    typedef StandardLocalAssemblerForIdentityOperatorOnSurface<
        BasisValueType, ResultType, GeometryFactory>
        LocalAssemblerForIdentityOperator_;
    typedef StandardLocalAssemblerForGridFunctionsOnSurfaces<
        BasisValueType, ResultType, GeometryFactory>
        LocalAssemblerForGridFunctions_;
    typedef StandardEvaluatorForIntegralOperators<
        BasisValueType, ResultType, GeometryFactory>
        EvaluatorForIntegralOperators_;

public:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& testBases,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Expression<BasisValueType>& testExpression,
            const Kernel<ResultType>& kernel,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const {
        return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                    new LocalAssemblerForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, kernel, trialExpression,
                        openClHandler, cacheSingularIntegrals,
                        m_accuracyOptions));
    }

    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const {
        throw std::runtime_error("StandardLocalAssemblerFactoryForOperatorsOnSurfaces::"
                                 "make(): collocation mode not implemented yet.");
    }

    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& testBases,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Expression<BasisValueType>& testExpression,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        return std::auto_ptr<LocalAssemblerForOperators<ResultType> >(
                    new LocalAssemblerForIdentityOperator_(
                        geometryFactory, rawGeometry,
                        testBases, trialBases,
                        testExpression, trialExpression,
                        openClHandler));
    }

    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
    throw std::runtime_error("StandardLocalAssemblerFactoryForOperatorsOnSurfaces::"
                             "make(): collocation mode not implemented yet.");
    }

    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& testBases,
            const Expression<BasisValueType>& testExpression,
            const Function<ResultType>& function,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        return std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >(
                    new LocalAssemblerForGridFunctions_(
                        geometryFactory, rawGeometry,
                        testBases,
                        testExpression, function,
                        openClHandler));
    }

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<BasisValueType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        return std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >(
                    new EvaluatorForIntegralOperators_(
                        geometryFactory, rawGeometry,
                        trialBases,
                        kernel, trialExpression, argumentLocalCoefficients,
                        openClHandler,
                        m_accuracyOptions.singleRegular));
    }

private:
    AccuracyOptions m_accuracyOptions;
};

} // namespace Fiber

#endif
