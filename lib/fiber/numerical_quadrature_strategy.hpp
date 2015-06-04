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

namespace Fiber {

template <typename BasisFunctionType> class QuadratureDescriptorSelectorFactory;
template <typename CoordinateType> class DoubleQuadratureRuleFamily;
template <typename CoordinateType> class SingleQuadratureRuleFamily;

/** \ingroup quadrature
 *  \brief Base class for NumericalQuadratureStrategy.
 *
 * This is the base class of the default quadrature strategy available
 * in BEM++. This quadrature strategy evaluates integrals by numerical
 * quadrature.
 *
 * The process of selecting a quadrature rule for the evaluation of a
 * particular integral can be customized at different levels of
 * generality. The choice of quadrature rule is done in two steps:
 *
 * 1. A *quadrature descriptor selector* is given the index of the
 *    element, or the pair of elements, over which integration should
 *    be done. It determines the desired order of accuracy of the
 *    quadrature rule and, for integrals over pairs of elements, the
 *    configuration of the two elements, i.e. whether they are
 *    coincident, adjacent or disjoint. These pieces of information
 *    are stored in a *quadrature descriptor*.
 *
 * 2. A *quadrature rule family* is given a quadrature descriptor and
 *    determines the points and weights of the quadrature rule to be
 *    applied.
 *
 * By default, NumericalQuadratureStrategy uses quadrature descriptor
 * selectors being instances of the classes
 * DefaultQuadratureDescriptorSelectorForIntegralOperators,
 * DefaultQuadratureDescriptorSelectorForLocalOperators and
 * DefaultQuadratureDescriptorSelectorForGridFunctions. You can make
 * it use different selectors by passing a custom
 * QuadratureDescriptorSelectorFactory object to the constructor of
 * NumericalQuadratureStrategy.
 *
 * The default quadrature descriptor selectors are customizable: you
 * can control the choice of quadrature orders using an
 * AccuracyOptionsEx options passed to another constructor of
 * NumericalQuadratureStrategy.
 *
 * By default, NumericalQuadratureStrategy uses the quadrature rule
 * families being instances of DefaultDoubleQuadratureRuleFamily and
 * DefaultSingleQuadratureRuleFamily. These use Gaussian quadrature
 * for regular integrals and the Sauter-Schwab quadrature rules for
 * singular integrals. If you wish, you can subclass
 * DoubleQuadratureRuleFamily and/or SingleQuadratureRuleFamily and
 * pass their instances to a NumericalQuadratureStrategy contructor.
 */
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable>
class NumericalQuadratureStrategyBase
    : public QuadratureStrategy<BasisFunctionType, ResultType, GeometryFactory,
                                Enable> {
public:
  typedef QuadratureStrategy<BasisFunctionType, ResultType, GeometryFactory,
                             Enable> Base;
  typedef typename Base::CoordinateType CoordinateType;

  NumericalQuadratureStrategyBase(
      const shared_ptr<const QuadratureDescriptorSelectorFactory<
          BasisFunctionType>> &quadratureDescriptorSelectorFactory,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> &
          singleQuadratureRuleFamily,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          doubleQuadratureRuleFamily);

public:
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
      const shared_ptr<const OpenClHandler> &openClHandler) const;

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
      const shared_ptr<const OpenClHandler> &openClHandler) const;

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
      const shared_ptr<const CollectionOfKernels<CoordinateType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, CoordinateType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals) const;

  virtual std::unique_ptr<LocalAssemblerForGridFunctions<ResultType>>
  makeAssemblerForGridFunctionsImplRealUserFunction(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const Function<CoordinateType>> &function,
      const shared_ptr<const OpenClHandler> &openClHandler) const;

  virtual std::unique_ptr<EvaluatorForIntegralOperators<ResultType>>
  makeEvaluatorForIntegralOperatorsImplRealKernel(
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
      const ParallelizationOptions &parallelizationOptions) const;

  virtual std::unique_ptr<LocalAssemblerForPotentialOperators<ResultType>>
  makeAssemblerForPotentialOperatorsImplRealKernel(
      const Matrix<CoordinateType> &evaluationPoints,
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
      VerbosityLevel::Level verbosityLevel) const;

protected:
  shared_ptr<const QuadratureDescriptorSelectorFactory<BasisFunctionType>>
  quadratureDescriptorSelectorFactory() const;
  shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
  singleQuadratureRuleFamily() const;
  shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>>
  doubleQuadratureRuleFamily() const;

private:
  shared_ptr<const QuadratureDescriptorSelectorFactory<BasisFunctionType>>
      m_quadratureDescriptorSelectorFactory;
  shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
      m_singleQuadratureRuleFamily;
  shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>>
      m_doubleQuadratureRuleFamily;
};

// Complex ResultType
/** \ingroup quadrature
 * \brief Quadrature strategy according to which integrals are
 * evaluated by numerical quadrature.
 *
 * A quadrature strategy provides functions constructing local assemblers used
 * to discretize boundary operators and user-defined functions. A particular
 * quadrature strategy determines how the integrals involved in this
 * discretization are evaluated.
 *
 * This is the default quadrature strategy available in BEM++. In this
 * quadrature strategy integrals are evaluated by numerical
 * quadrature.
 *
 * The process of selecting a quadrature rule for the evaluation of a
 * particular integral can be customized at different levels of
 * generality. The choice of quadrature rule is done in two steps:
 *
 * 1. A *quadrature descriptor selector* is given the index of the
 *    element, or the pair of elements, over which integration should
 *    be done. It determines the desired order of accuracy of the
 *    quadrature rule and, for integrals over pairs of elements, the
 *    configuration of the two elements, i.e. whether they are
 *    coincident, adjacent or disjoint. These pieces of information
 *    are stored in a *quadrature descriptor*.
 *
 * 2. A *quadrature rule family* is given a quadrature descriptor and
 *    determines the points and weights of the quadrature rule to be
 *    applied.
 *
 * By default, NumericalQuadratureStrategy uses quadrature descriptor
 * selectors being instances of the classes
 * DefaultQuadratureDescriptorSelectorForIntegralOperators,
 * DefaultQuadratureDescriptorSelectorForLocalOperators and
 * DefaultQuadratureDescriptorSelectorForGridFunctions. You can make
 * it use different selectors by passing a custom
 * QuadratureDescriptorSelectorFactory object to the constructor of
 * NumericalQuadratureStrategy.
 *
 * The default quadrature descriptor selectors are customizable: you
 * can control the choice of quadrature orders using an
 * AccuracyOptionsEx options passed to another constructor of
 * NumericalQuadratureStrategy.
 *
 * By default, NumericalQuadratureStrategy uses the quadrature rule
 * families being instances of DefaultDoubleQuadratureRuleFamily and
 * DefaultSingleQuadratureRuleFamily. These use Gaussian quadrature
 * for regular integrals and the Sauter-Schwab quadrature rules (*) for
 * singular integrals. If you wish, you can subclass
 * DoubleQuadratureRuleFamily and/or SingleQuadratureRuleFamily and
 * pass their instances to a NumericalQuadratureStrategy contructor.
 *
 * (*) S. Sauter, Ch. Schwab, "Boundary Element Methods" (2010).
 */
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory, typename Enable = void>
class NumericalQuadratureStrategy
    : public NumericalQuadratureStrategyBase<BasisFunctionType, ResultType,
                                             GeometryFactory, Enable> {
  typedef NumericalQuadratureStrategyBase<BasisFunctionType, ResultType,
                                          GeometryFactory, Enable> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use the default
   * quadrature descriptor selector factory and the default accuracy
   * options. */
  NumericalQuadratureStrategy();

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use the default
   * quadrature descriptor selector factory with custom accuracy
   * options and the default quadrature rule families. */
  explicit NumericalQuadratureStrategy(
      const AccuracyOptionsEx &accuracyOptions);

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use the default
   * quadrature descriptor selector factory with custom accuracy
   * options and custom quadrature rule families. */
  NumericalQuadratureStrategy(
      const AccuracyOptionsEx &accuracyOptions,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> &
          singleQuadratureRuleFamily,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          doubleQuadratureRuleFamily);

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use a custom
   * quadrature descriptor selector factory and quadrature rule families.
   */
  NumericalQuadratureStrategy(
      const shared_ptr<const QuadratureDescriptorSelectorFactory<
          BasisFunctionType>> &quadratureDescriptorSelectorFactory,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> &
          singleQuadratureRuleFamily,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          doubleQuadratureRuleFamily);

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
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals) const;

  virtual std::unique_ptr<LocalAssemblerForGridFunctions<ResultType>>
  makeAssemblerForGridFunctionsImplComplexUserFunction(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const Function<ResultType>> &function,
      const shared_ptr<const OpenClHandler> &openClHandler) const;

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
      const ParallelizationOptions &parallelizationOptions) const;

  virtual std::unique_ptr<LocalAssemblerForPotentialOperators<ResultType>>
  makeAssemblerForPotentialOperatorsImplComplexKernel(
      const Matrix<CoordinateType> &evaluationPoints,
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
      VerbosityLevel::Level verbosityLevel) const;
};

/** \cond ENABLE_IFS */
// RealResultType
template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
class NumericalQuadratureStrategy<
    BasisFunctionType, ResultType, GeometryFactory,
    typename boost::enable_if<boost::is_same<
        ResultType, typename ScalarTraits<ResultType>::RealType>>::type>
    : public NumericalQuadratureStrategyBase<
          BasisFunctionType, ResultType, GeometryFactory,
          typename boost::enable_if<boost::is_same<
              ResultType, typename ScalarTraits<ResultType>::RealType>>::type> {
  typedef typename boost::enable_if<boost::is_same<
      ResultType, typename ScalarTraits<ResultType>::RealType>>::type Enable;
  typedef NumericalQuadratureStrategyBase<BasisFunctionType, ResultType,
                                          GeometryFactory, Enable> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use the default
   * quadrature descriptor selector factory and the default accuracy
   * options. */
  NumericalQuadratureStrategy();

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use the default
   * quadrature descriptor selector factory with custom accuracy
   * options and the default quadrature rule families. */
  explicit NumericalQuadratureStrategy(
      const AccuracyOptionsEx &accuracyOptions);

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use a custom
   * quadrature descriptor selector factory and quadrature rule families.
   */
  explicit NumericalQuadratureStrategy(
      const shared_ptr<const QuadratureDescriptorSelectorFactory<
          BasisFunctionType>> &quadratureDescriptorSelectorFactory,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> &
          singleQuadratureRuleFamily,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          doubleQuadratureRuleFamily);

  /** \brief Construct a numerical quadrature strategy.
   *
   * This constructor makes the newly created object use the default
   * quadrature descriptor selector factory with custom accuracy
   * options and custom quadrature rule families. */
  explicit NumericalQuadratureStrategy(
      const AccuracyOptionsEx &accuracyOptions,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> &
          singleQuadratureRuleFamily,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          doubleQuadratureRuleFamily);
};
/** \endcond */

} // namespace Fiber

#include "numerical_quadrature_strategy_imp.hpp"

#endif
