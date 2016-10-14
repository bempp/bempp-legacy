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

#ifndef fiber_cuda_separable_numerical_test_kernel_trial_integrator_hpp
#define fiber_cuda_separable_numerical_test_kernel_trial_integrator_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegrator;
template <typename BasisFunctionType> class Shapeset;
/** \endcond */

/** \brief Integration over pairs of elements on tensor-product point grids
 * with CUDA. */
template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class CudaSeparableNumericalTestKernelTrialIntegrator
    : public TestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                       ResultType> {
public:

  typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
      Base;
  typedef typename Base::CoordinateType CoordinateType;

  CudaSeparableNumericalTestKernelTrialIntegrator(
      const Matrix<CoordinateType> &localTestQuadPoints,
      const Matrix<CoordinateType> &localTrialQuadPoints,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights);

  virtual ~CudaSeparableNumericalTestKernelTrialIntegrator();

  virtual void integrate(CallVariant callVariant,
                         const std::vector<int> &elementIndicesA,
                         int elementIndexB,
                         const Shapeset<BasisFunctionType> &basisA,
                         const Shapeset<BasisFunctionType> &basisB,
                         LocalDofIndex localDofIndexB,
                         const std::vector<Matrix<ResultType> *> &result) const;
private:

  /** brief Evaluate the kernels on a grid of test and trial points
   * \param[in] testPoints Global coordinates of points on the test element
   * \param[in] trialPoints Global coordinates of point on the trial element
   * \param[out] kernelValues Results
   */
//  void evaluateKernel(thrust::device_vector<double> &testPoints,
//                      thrust::device_vector<double> &trialPoints,
//                      thrust::device_vector<double> &kernelValues);
//
//  void evaluateIntegral(const unsigned int testElemCount,
//                        const unsigned int trialElemCount,
//                        thrust::device_vector<double> &kernelValues,
//                        thrust::device_vector<double> &testBasisData,
//                        thrust::device_vector<double> &trialBasisData,
//                        thrust::device_vector<double> &testQuadWeights,
//                        thrust::device_vector<double> &trialQuadWeights,
//                        thrust::device_vector<double> &result);

  Matrix<CoordinateType> m_localTestQuadPoints;
  Matrix<CoordinateType> m_localTrialQuadPoints;
  std::vector<CoordinateType> m_testQuadWeights;
  std::vector<CoordinateType> m_trialQuadWeights;
};

} // namespace Fiber

#endif
