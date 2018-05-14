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

#ifndef fiber_cuda_evaluate_laplace_3d_adjoint_double_layer_potential_integral_cuh
#define fiber_cuda_evaluate_laplace_3d_adjoint_double_layer_potential_integral_cuh

#include "cuda.hpp"

#include "../common/scalar_traits.hpp"

#include <device_launch_parameters.h>

namespace Fiber {

template<typename ValueType>
__device__ void Laplace3dAdjointDoubleLayerPotentialKernel(
    typename ScalarTraits<ValueType>::RealType testPointCoo[3],
    typename ScalarTraits<ValueType>::RealType trialPointCoo[3],
    typename ScalarTraits<ValueType>::RealType testElemNormal[3],
    ValueType &result) {

  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  CoordinateType numeratorSum = 0., distanceSq = 0.;
#pragma unroll
  for (int coordIndex = 0; coordIndex < 3; ++coordIndex) {
    CoordinateType diff =
        testPointCoo[coordIndex] - trialPointCoo[coordIndex];
    distanceSq += diff * diff;
    numeratorSum += diff * testElemNormal[coordIndex];
  }
  CoordinateType distance = sqrt(distanceSq);
  result = -numeratorSum / (static_cast<CoordinateType>(4. * M_PI) *
                            distanceSq * distance);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
CudaEvaluateLaplace3dAdjointDoubleLayerPotentialIntegralFunctorCached(
    const int elemPairIndexBegin, const unsigned int elemPairCount,
    const unsigned int testIndexCount,
    const int* __restrict__ testIndices, const int* __restrict__ trialIndices,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const unsigned int testDofCount,
    const BasisFunctionType* __restrict__ testBasisValues,
    const unsigned int trialDofCount,
    const BasisFunctionType* __restrict__ trialBasisValues,
    const unsigned int testElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* testElemNormals,
    const typename ScalarTraits<BasisFunctionType>::RealType* testIntegrationElements,
    const unsigned int trialElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialGlobalPoints,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialIntegrationElements,
    ResultType* __restrict__ result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < elemPairCount) {

    const int coordCount = 3;

    // Determine trial and test element indices
    const int elemPairIndex = elemPairIndexBegin + i;
    const int trialElemPosition = trialIndices[elemPairIndex / testIndexCount];
    const int testElemPosition = testIndices[elemPairIndex % testIndexCount];

    // Gather integration elements
    const CoordinateType trialIntegrationElement =
        trialIntegrationElements[trialElemPosition];
    const CoordinateType testIntegrationElement =
        testIntegrationElements[testElemPosition];

    // Get normals
    CoordinateType testElemNormal[coordCount];
#pragma unroll
    for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      testElemNormal[coordIndex] =
          testElemNormals[coordIndex * testElemCount + testElemPosition];
    }

    // Evaluate kernel
    KernelType kernelValues[6 * 6];
    CoordinateType trialPointCoo[coordCount], testPointCoo[coordCount];
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
        KernelType kernelValue;
#pragma unroll
        for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
          trialPointCoo[coordIndex] =
              trialGlobalPoints[coordIndex * trialPointCount * trialElemCount
                            + trialPoint * trialElemCount
                            + trialElemPosition];
          testPointCoo[coordIndex] =
              testGlobalPoints[coordIndex * testPointCount * testElemCount
                           + testPoint * testElemCount
                           + testElemPosition];
        }
        Laplace3dAdjointDoubleLayerPotentialKernel(
            testPointCoo, trialPointCoo, testElemNormal,
            kernelValue);
        kernelValues[trialPoint * testPointCount + testPoint] = kernelValue;
      }
    }

    // Perform numerical integration
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        ResultType sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * constTrialQuadWeights[2][trialPoint];
          ResultType partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * constTestQuadWeights[2][testPoint];
            partialSum +=
                kernelValues[trialPoint * testPointCount + testPoint]
                * testBasisValues[testDof * testPointCount + testPoint]
                * trialBasisValues[trialDof * trialPointCount + trialPoint]
                * testWeight;
          }
          sum += partialSum * trialWeight;
        }
        const size_t index = trialDof * testDofCount * elemPairCount
                             + testDof * elemPairCount
                             + i;
        result[index] = sum;
      }
    }
  }
}

} // namespace Fiber

#endif
