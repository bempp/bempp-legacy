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

#ifndef fiber_default_evaluator_for_integral_operators_hpp
#define fiber_default_evaluator_for_integral_operators_hpp

#include "../common/common.hpp"

#include "evaluator_for_integral_operators.hpp"

#include "collection_of_2d_arrays.hpp"
#include "parallelization_options.hpp"
#include "quadrature_options.hpp"
#include "types.hpp"
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
class QuadratureOptions;
template <typename ValueType> class Shapeset;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
template <typename CoordinateType> class RawGridGeometry;
class OpenClHandler;
template <typename BasisFunctionType>
class QuadratureDescriptorSelectorForPotentialOperators;
template <typename CoordinateType> class SingleQuadratureRuleFamily;
/** \endcond */

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class DefaultEvaluatorForIntegralOperators
    : public EvaluatorForIntegralOperators<ResultType> {
public:
  typedef EvaluatorForIntegralOperators<ResultType> Base;
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::Region Region;

  DefaultEvaluatorForIntegralOperators(
      const shared_ptr<const GeometryFactory> &geometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfKernels<KernelType>> &kernels,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const KernelTrialIntegral<BasisFunctionType, KernelType,
                                                 ResultType>> &integral,
      const shared_ptr<const std::vector<std::vector<ResultType>>> &
          argumentLocalCoefficients,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      const shared_ptr<const QuadratureDescriptorSelectorForPotentialOperators<
          BasisFunctionType>> &quadDescSelector,
      const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>> &
          quadRuleFamily);

  virtual void evaluate(Region region, const Matrix<CoordinateType> &points,
                        Matrix<ResultType> &result) const;

private:
  void cacheTrialData();
  void calcTrialData(Region region, int kernelTrialGeomDeps,
                     GeometricalData<CoordinateType> &trialGeomData,
                     CollectionOf2dArrays<ResultType> &trialExprValues,
                     std::vector<CoordinateType> &weights) const;

private:
  const shared_ptr<const GeometryFactory> m_geometryFactory;
  const shared_ptr<const RawGridGeometry<CoordinateType>> m_rawGeometry;
  const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
      m_trialShapesets;
  const shared_ptr<const CollectionOfKernels<KernelType>> m_kernels;
  const shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
      m_trialTransformations;
  const shared_ptr<const KernelTrialIntegral<BasisFunctionType, KernelType,
                                             ResultType>> m_integral;
  const shared_ptr<const std::vector<std::vector<ResultType>>>
      m_argumentLocalCoefficients;
  const shared_ptr<const OpenClHandler> m_openClHandler;
  const ParallelizationOptions m_parallelizationOptions;
  const shared_ptr<const QuadratureDescriptorSelectorForPotentialOperators<
      BasisFunctionType>> m_quadDescSelector;
  const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
      m_quadRuleFamily;

  Fiber::GeometricalData<CoordinateType> m_nearFieldTrialGeomData;
  Fiber::GeometricalData<CoordinateType> m_farFieldTrialGeomData;
  CollectionOf2dArrays<ResultType> m_nearFieldTrialTransfValues;
  CollectionOf2dArrays<ResultType> m_farFieldTrialTransfValues;
  std::vector<CoordinateType> m_nearFieldWeights;
  std::vector<CoordinateType> m_farFieldWeights;
};

} // namespace Fiber

#include "default_evaluator_for_integral_operators_imp.hpp"

#endif
