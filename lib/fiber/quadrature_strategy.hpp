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

#ifndef fiber_quadrature_strategy_hpp
#define fiber_quadrature_strategy_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"
#include "shared_ptr.hpp"
#include "verbosity_level.hpp"

#include "../common/armadillo_fwd.hpp"
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <memory>

namespace Fiber {

/** \cond FORWARD_DECL */
class ParallelizationOptions;
class OpenClHandler;

template <typename ValueType> class Shapeset;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class Function;
template <typename ValueType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
template <typename BasisFunctionType, typename ResultType>
class TestTrialIntegral;
template <typename CoordinateType> class RawGridGeometry;

template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename ResultType> class LocalAssemblerForLocalOperators;
template <typename ResultType> class LocalAssemblerForPotentialOperators;
template <typename ResultType> class LocalAssemblerForGridFunctions;
template <typename ResultType> class EvaluatorForIntegralOperators;
/** \endcond */

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class QuadratureStrategyBase {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  virtual ~QuadratureStrategyBase() {}

  /** \brief Allocate a Galerkin-mode local assembler for an integral operator
      with real kernel. */
  std::unique_ptr<LocalAssemblerForIntegralOperators<ResultType>>
  makeAssemblerForIntegralOperators(
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals) const {
    return this->makeAssemblerForIntegralOperatorsImplRealKernel(
        testGeometryFactory, trialGeometryFactory, testRawGeometry,
        trialRawGeometry, testShapesets, trialShapesets, testTransformations,
        kernels, trialTransformations, integral, openClHandler,
        parallelizationOptions, verbosityLevel, cacheSingularIntegrals);
  }

  /** \brief Allocate a Galerkin-mode local assembler for the identity operator.
   *
   *  \deprecated This method is deprecated. Use the more general
   *  makeAssemblerForLocalOperators() method, passing an appropriate
   *  TestTrialIntegral object. */
  virtual std::unique_ptr<LocalAssemblerForLocalOperators<ResultType>>
  makeAssemblerForIdentityOperators(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const OpenClHandler> &openClHandler) const = 0;

  /** \brief Allocate a Galerkin-mode local assembler for a local operator. */
  virtual std::unique_ptr<LocalAssemblerForLocalOperators<ResultType>>
  makeAssemblerForLocalOperators(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestTrialIntegral<BasisFunctionType, ResultType>> &
          integral,
      const shared_ptr<const OpenClHandler> &openClHandler) const = 0;

  /** \brief Allocate a local assembler for calculations of the projections
    of functions from a given space on a Fiber::Function. */
  std::unique_ptr<LocalAssemblerForGridFunctions<ResultType>>
  makeAssemblerForGridFunctions(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const Function<ResultType>> &function,
      const shared_ptr<const OpenClHandler> &openClHandler) const {
    return this->makeAssemblerForGridFunctionsImplRealUserFunction(
        geometryFactory, rawGeometry, testShapesets, testTransformations,
        function, openClHandler);
  }

  /** \brief Allocate an evaluator for an integral operator with real kernel
    applied to a grid function. */
  std::unique_ptr<EvaluatorForIntegralOperators<ResultType>>
  makeEvaluatorForIntegralOperators(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const std::vector<std::vector<ResultType>>> &
          argumentLocalCoefficients,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions) const {
    return this->makeEvaluatorForIntegralOperatorsImplRealKernel(
        geometryFactory, rawGeometry, trialShapesets, kernels,
        trialTransformations, integral, argumentLocalCoefficients,
        openClHandler, parallelizationOptions);
  }

  /** \brief Allocate a local assembler for a potential operator with real
   * kernel. */
  std::unique_ptr<LocalAssemblerForPotentialOperators<ResultType>>
  makeAssemblerForPotentialOperators(
      const arma::Mat<CoordinateType> &evaluationPoints,
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel) const {
    return this->makeAssemblerForPotentialOperatorsImplRealKernel(
        evaluationPoints, geometryFactory, rawGeometry, trialShapesets, kernels,
        trialTransformations, integral, openClHandler, parallelizationOptions,
        verbosityLevel);
  }

private:
  virtual std::unique_ptr<LocalAssemblerForIntegralOperators<ResultType>>
  makeAssemblerForIntegralOperatorsImplRealKernel(
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernel,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel,
      bool cacheSingularIntegrals) const = 0;

  virtual std::unique_ptr<LocalAssemblerForGridFunctions<ResultType>>
  makeAssemblerForGridFunctionsImplRealUserFunction(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const Function<CoordinateType>> &function,
      const shared_ptr<const OpenClHandler> &openClHandler) const = 0;

  virtual std::unique_ptr<EvaluatorForIntegralOperators<ResultType>>
  makeEvaluatorForIntegralOperatorsImplRealKernel(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernel,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const std::vector<std::vector<ResultType>>> &
          argumentLocalCoefficients,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions) const = 0;

  virtual std::unique_ptr<LocalAssemblerForPotentialOperators<ResultType>>
  makeAssemblerForPotentialOperatorsImplRealKernel(
      const arma::Mat<CoordinateType> &evaluationPoints,
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel) const = 0;
};

/** \ingroup quadrature
 *  \brief Base class for quadrature strategies.
 *
 *  A quadrature strategy provides functions constructing local assemblers used
 *  to discretize boundary operators and user-defined functions. A particular
 *  quadrature strategy determines how the integrals involved in this
 *  discretization are evaluated. */
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable = void>
class QuadratureStrategy : public QuadratureStrategyBase<
                               BasisFunctionType, ResultType, GeometryFactory> {
  typedef QuadratureStrategyBase<BasisFunctionType, ResultType, GeometryFactory>
  Base;

public:
  typedef typename Base::CoordinateType CoordinateType;

  using Base::makeAssemblerForGridFunctions;
  using Base::makeAssemblerForIntegralOperators;
  using Base::makeAssemblerForPotentialOperators;
  using Base::makeEvaluatorForIntegralOperators;

  /** \brief Allocate a Galerkin-mode local assembler for an integral operator
      with complex kernel. */
  std::unique_ptr<LocalAssemblerForIntegralOperators<ResultType>>
  makeAssemblerForIntegralOperators(
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfKernels<ResultType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, ResultType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals) const {
    return this->makeAssemblerForIntegralOperatorsImplComplexKernel(
        testGeometryFactory, trialGeometryFactory, testRawGeometry,
        trialRawGeometry, testShapesets, trialShapesets, testTransformations,
        kernels, trialTransformations, integral, openClHandler,
        parallelizationOptions, verbosityLevel, cacheSingularIntegrals);
  }

  /** \brief Allocate a local assembler for calculations of the projections
    of complex-valued functions from a given space on a Fiber::Function. */
  std::unique_ptr<LocalAssemblerForGridFunctions<ResultType>>
  makeAssemblerForGridFunctions(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const Function<ResultType>> &function,
      const shared_ptr<const OpenClHandler> &openClHandler) const {
    return this->makeAssemblerForGridFunctionsImplComplexUserFunction(
        geometryFactory, rawGeometry, testShapesets, testTransformations,
        function, openClHandler);
  }

  /** \brief Allocate an evaluator for an integral operator with
    complex-valued kernel applied to a grid function. */
  std::unique_ptr<EvaluatorForIntegralOperators<ResultType>>
  makeEvaluatorForIntegralOperators(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<ResultType>> &kernel,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType,
                                                 ResultType>> &integral,
      const shared_ptr<const std::vector<std::vector<ResultType>>> &
          argumentLocalCoefficients,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions) const {
    return this->makeEvaluatorForIntegralOperatorsImplComplexKernel(
        geometryFactory, rawGeometry, trialShapesets, kernel,
        trialTransformations, integral, argumentLocalCoefficients,
        openClHandler, parallelizationOptions);
  }

  /** \brief Allocate a local assembler for a potential operator with complex
   *  kernel. */
  std::unique_ptr<LocalAssemblerForPotentialOperators<ResultType>>
  makeAssemblerForPotentialOperators(
      const arma::Mat<CoordinateType> &evaluationPoints,
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<ResultType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType,
                                                 ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel) const {
    return this->makeAssemblerForPotentialOperatorsImplComplexKernel(
        evaluationPoints, geometryFactory, rawGeometry, trialShapesets, kernels,
        trialTransformations, integral, openClHandler, parallelizationOptions,
        verbosityLevel);
  }

private:
  virtual std::unique_ptr<LocalAssemblerForIntegralOperators<ResultType>>
  makeAssemblerForIntegralOperatorsImplComplexKernel(
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfKernels<ResultType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, ResultType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel,
      bool cacheSingularIntegrals) const = 0;

  virtual std::unique_ptr<LocalAssemblerForGridFunctions<ResultType>>
  makeAssemblerForGridFunctionsImplComplexUserFunction(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const Function<ResultType>> &function,
      const shared_ptr<const OpenClHandler> &openClHandler) const = 0;

  virtual std::unique_ptr<EvaluatorForIntegralOperators<ResultType>>
  makeEvaluatorForIntegralOperatorsImplComplexKernel(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<ResultType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType,
                                                 ResultType>> &integral,
      const shared_ptr<const std::vector<std::vector<ResultType>>> &
          argumentLocalCoefficients,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions) const = 0;

  virtual std::unique_ptr<LocalAssemblerForPotentialOperators<ResultType>>
  makeAssemblerForPotentialOperatorsImplComplexKernel(
      const arma::Mat<CoordinateType> &evaluationPoints,
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<ResultType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<BasisFunctionType, ResultType,
                                                 ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel) const = 0;
};

/** \cond ENABLE_IFS */
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class QuadratureStrategy<
    BasisFunctionType, ResultType, GeometryFactory,
    typename boost::enable_if<boost::is_same<
        ResultType, typename ScalarTraits<ResultType>::RealType>>::
        type> : public QuadratureStrategyBase<BasisFunctionType, ResultType,
                                              GeometryFactory> {
  typedef QuadratureStrategyBase<BasisFunctionType, ResultType, GeometryFactory>
  Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
};
/** \endcond  */

} // namespace Fiber

#endif
