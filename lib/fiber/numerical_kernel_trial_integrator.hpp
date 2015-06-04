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

#ifndef fiber_numerical_kernel_trial_integrator_hpp
#define fiber_numerical_kernel_trial_integrator_hpp

#include "../common/common.hpp"

#include "kernel_trial_integrator.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename CoordinateType> class RawGridGeometry;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
/** \endcond */

/** \brief Integration over pairs of elements on tensor-product point grids. */
template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class NumericalKernelTrialIntegrator
    : public KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> {
public:
  typedef KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> Base;
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::PointElementIndexPair PointElementIndexPair;

  NumericalKernelTrialIntegrator(
      const Matrix<CoordinateType> &localQuadPoints,
      const std::vector<CoordinateType> quadWeights,
      const Matrix<CoordinateType> &points,
      const GeometryFactory &geometryFactory,
      const RawGridGeometry<CoordinateType> &rawGeometry,
      const CollectionOfKernels<KernelType> &kernels,
      const CollectionOfShapesetTransformations<CoordinateType>
          &trialTransformations,
      const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
          &integral);

  virtual void integrate(const std::vector<int> &pointIndices,
                         int trialElementIndex,
                         const Shapeset<BasisFunctionType> &trialShapeset,
                         LocalDofIndex localTrialDofIndex,
                         const std::vector<Matrix<ResultType> *> &result) const;

  virtual void integrate(int pointIndex, int componentIndex,
                         const std::vector<int> &trialElementIndices,
                         const Shapeset<BasisFunctionType> &trialShapeset,
                         const std::vector<Matrix<ResultType> *> &result) const;

  virtual void
  integrate(const std::vector<PointElementIndexPair> &pointElementIndexPairs,
            const Shapeset<BasisFunctionType> &trialShapeset,
            const std::vector<Matrix<ResultType> *> &result) const;

private:
  /** \cond PRIVATE */
  Matrix<CoordinateType> m_localQuadPoints;
  std::vector<CoordinateType> m_quadWeights;

  const Matrix<CoordinateType> &m_points;

  const GeometryFactory &m_geometryFactory;
  const RawGridGeometry<CoordinateType> &m_rawGeometry;

  const CollectionOfKernels<KernelType> &m_kernels;
  const CollectionOfShapesetTransformations<CoordinateType>
      &m_trialTransformations;
  const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
      &m_integral;
  /** \endcond */
};

} // namespace Fiber

#include "numerical_kernel_trial_integrator_imp.hpp"

#endif
