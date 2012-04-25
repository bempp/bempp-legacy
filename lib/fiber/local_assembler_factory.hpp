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

template <typename BasisValueType, typename ResultType,
          typename GeometryFactory>
class LocalAssemblerFactory
{
public:
    // Type of the output of integral operators. Note that identity operators
    // output numbers of type BasisValueType.
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    virtual ~LocalAssemblerFactory() {}

    /** @name Local assemblers for integral operators
        @{ */

    // TODO (important): overload for Kernel<ResultType> and
    // Kernel<ScalarTraits<ResultType> >::RealType
    /** \brief Allocate a Galerkin-mode local assembler for an integral operator. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& testBases,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Expression<BasisValueType>& testExpression,
            const Kernel<ResultType>& kernel,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** \brief Allocate a collocation-mode local assembler for an integral operator.

        Used also for evaluation of the integral operator at arbitrary points. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler,
            bool cacheSingularIntegrals) const = 0;

    /** @}
        @name Local assemblers for the identity operator
        @{ */

    /** \brief Allocate a Galerkin-mode local assembler for the identity operator. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& testBases,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Expression<BasisValueType>& testExpression,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;

    /** \brief Allocate a collocation-mode local assembler for an identity operator.

        Used also for evaluation of the identity operator at arbitrary points. */
    virtual std::auto_ptr<LocalAssemblerForOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Expression<BasisValueType>& trialExpression,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;

    /** @}
        @name Local assemblers for grid functions
        @{ */

    /** \brief Allocate a local assembler for calculations of the projections
      of functions from a given space on a Fiber::Function. */
    virtual std::auto_ptr<LocalAssemblerForGridFunctions<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& testBases,
            const Expression<BasisValueType>& testExpression,
            const Function<ResultType>& function,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;

    /** @}
        @name Evaluators for integral operators
        @{ */

    /** \brief Allocate an evaluator for an integral operator applied to a
      grid function. */
    virtual std::auto_ptr<EvaluatorForIntegralOperators<ResultType> > make(
            const GeometryFactory& geometryFactory,
            const RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Basis<BasisValueType>*>& trialBases,
            const Kernel<ResultType>& kernel,
            const Expression<BasisValueType>& trialExpression,
            const std::vector<std::vector<ResultType> >& argumentLocalCoefficients,
            const OpenClHandler<CoordinateType, int>& openClHandler) const = 0;

    /** @} */
};

} // namespace Fiber

#endif
