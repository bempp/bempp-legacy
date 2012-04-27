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

#ifndef fiber_integration_manager_factory_hpp
#define fiber_integration_manager_factory_hpp

#include "scalar_traits.hpp"

#include <armadillo>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <memory>

namespace Fiber
{
template <typename ValueType, typename IndexType> class OpenClHandler;

template <typename ValueType> class Basis;
template <typename ValueType> class Expression;
template <typename ValueType> class Function;
template <typename ValueType> class Kernel;
template <typename CoordinateType> class RawGridGeometry;

template <typename ResultType> class LocalAssemblerForOperators;
template <typename ResultType> class LocalAssemblerForGridFunctions;
template <typename ResultType> class EvaluatorForIntegralOperators;

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class LocalAssemblerFactoryBase
{
public:
    // Type of the output of integral operators. Note that identity operators
    // output numbers of type BasisFunctionType.
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
private:
//    typedef mpl::set<ResultType, CoordinateType> AcceptedTypes;

public:
    virtual ~LocalAssemblerFactoryBase() {}

    /** \brief Allocate a Galerkin-mode local assembler for an integral operator
        with real kernel. */
    std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Kernel<CoordinateType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const {
        return this->makeAssemblerForIntegralOperatorsImplRealKernel(
                    geometryFactory, rawGeometry, testBases, trialBases,
                    testExpression, kernel, trialExpression, openClHandler,
                    cacheSingularIntegrals);
    }

    /** \brief Allocate a Galerkin-mode local assembler for the identity operator. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIdentityOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;

    /** \brief Allocate a local assembler for calculations of the projections
      of functions from a given space on a Fiber::Function. */
    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> >
    makeAssemblerForGridFunctions(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const Expression<CoordinateType>& testExpression,
            const Function<ResultType>& function,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;

    /** \brief Allocate an evaluator for an integral operator with real kernel
      applied to a grid function. */
    std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<CoordinateType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        return this->makeEvaluatorForIntegralOperatorsImplRealKernel(
                    geometryFactory, rawGeometry, trialBases,
                    kernel, trialExpression, argumentLocalCoefficients,
                    openClHandler);
    }

private:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperatorsImplRealKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Kernel<CoordinateType>& kernel, // !
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperatorsImplRealKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<CoordinateType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;
};

// complex ResultType
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable = void>
class LocalAssemblerFactory :
        public LocalAssemblerFactoryBase<BasisFunctionType, ResultType, GeometryFactory>
{
    typedef LocalAssemblerFactoryBase<BasisFunctionType, ResultType, GeometryFactory> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;

    using Base::makeAssemblerForIntegralOperators;
    using Base::makeEvaluatorForIntegralOperators;

    std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Kernel<ResultType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const {
        return this->makeAssemblerForIntegralOperatorsImplComplexKernel(
                    geometryFactory, rawGeometry, testBases, trialBases,
                    testExpression, kernel, trialExpression, openClHandler,
                    cacheSingularIntegrals);
    }

    /** \brief Allocate an evaluator for an integral operator with real kernel
      applied to a grid function. */
    std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperators(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const {
        return this->makeEvaluatorForIntegralOperatorsImplComplexKernel(
                    geometryFactory, rawGeometry, trialBases,
                    kernel, trialExpression, argumentLocalCoefficients,
                    openClHandler);
    }

private:
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> >
    makeAssemblerForIntegralOperatorsImplComplexKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Expression<CoordinateType>& testExpression,
            const Kernel<ResultType>& kernel, // !
            const Expression<CoordinateType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> >
    makeEvaluatorForIntegralOperatorsImplComplexKernel(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisFunctionType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<CoordinateType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;
};

// real ResultType
template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
class LocalAssemblerFactory<BasisFunctionType, ResultType, GeometryFactory,
    typename boost::enable_if<boost::is_same<ResultType, typename ScalarTraits<ResultType>::RealType> >::type > :
        public LocalAssemblerFactoryBase<BasisFunctionType, ResultType, GeometryFactory>
{
    typedef LocalAssemblerFactoryBase<BasisFunctionType, ResultType, GeometryFactory> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
};

} // namespace Fiber

#endif
