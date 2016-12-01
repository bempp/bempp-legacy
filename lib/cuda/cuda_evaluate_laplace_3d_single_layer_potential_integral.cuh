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

#ifndef fiber_cuda_evaluate_laplace_3d_single_layer_potential_integral_cuh
#define fiber_cuda_evaluate_laplace_3d_single_layer_potential_integral_cuh

#include "cuda.cuh"

#include "../common/scalar_traits.hpp"

#include <device_launch_parameters.h>

namespace Fiber {

template<typename ValueType>
__device__ void Laplace3dSingleLayerPotentialKernel(
    typename ScalarTraits<ValueType>::RealType testPointCoo[3],
    typename ScalarTraits<ValueType>::RealType trialPointCoo[3],
    ValueType &result) {

  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  ValueType sum = 0;
#pragma unroll
  for (int coordIndex = 0; coordIndex < 3; ++coordIndex) {
    ValueType diff = testPointCoo[coordIndex] - trialPointCoo[coordIndex];
    sum += diff * diff;
  }
  result = static_cast<CoordinateType>(1. / (4. * M_PI)) / sqrt(sum);

//  result = static_cast<CoordinateType>(1. / (4. * M_PI)) * rnorm3df(
//      testPointCoo[0] - trialPointCoo[0],
//      testPointCoo[1] - trialPointCoo[1],
//      testPointCoo[2] - trialPointCoo[2]);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
CudaEvaluateLaplace3dSingleLayerPotentialKernelFunctorCached(
    const int elemPairIndexBegin, const unsigned int elemPairCount,
    const unsigned int testIndexCount,
    const int* __restrict__ testIndices, const int* __restrict__ trialIndices,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const unsigned int testElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testGlobalPoints,
    const unsigned int trialElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialGlobalPoints,
    KernelType* __restrict__ kernelValues) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  // Each thread is working on one element pair
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < elemPairCount) {

    const int coordCount = 3;

    // Determine trial and test element indices
    const int elemPairIndex = elemPairIndexBegin + i;
    const int trialElemPosition = trialIndices[elemPairIndex / testIndexCount];
    const int testElemPosition = testIndices[elemPairIndex % testIndexCount];

    // Evaluate kernel
    KernelType kernelValue;
    CoordinateType trialPointCoo[coordCount], testPointCoo[coordCount];
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
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
        Laplace3dSingleLayerPotentialKernel(
            testPointCoo, trialPointCoo,
            kernelValue);
        const size_t index = trialPoint * testPointCount * elemPairCount
                             + testPoint * elemPairCount
                             + i;
        kernelValues[index] = kernelValue;
      }
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorKernelCached(
    const int elemPairIndexBegin, const unsigned int elemPairCount,
    const unsigned int testIndexCount,
    const int* __restrict__ testIndices, const int* __restrict__ trialIndices,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const unsigned int testDofCount,
    const BasisFunctionType* __restrict__ testBasisValues,
    const unsigned int trialDofCount,
    const BasisFunctionType* __restrict__ trialBasisValues,
    const unsigned int testElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testIntegrationElements,
    const unsigned int trialElemCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialIntegrationElements,
    const KernelType* __restrict__ kernelValues,
    ResultType* __restrict__ result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  // Each thread is working on one dof pair
  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < elemPairCount * trialDofCount * testDofCount) {

    // Determine trial and test element indices
    const int localElemPairIndex = i % elemPairCount;
    const int elemPairIndex = elemPairIndexBegin + localElemPairIndex;
    const int trialElemPosition = trialIndices[elemPairIndex / testIndexCount];
    const int testElemPosition = testIndices[elemPairIndex % testIndexCount];
    const int trialDof = i / (elemPairCount * testDofCount);
    const int testDof = (i % (elemPairCount * testDofCount)) / elemPairCount;

    // Gather integration elements
    const CoordinateType trialIntegrationElement =
        trialIntegrationElements[trialElemPosition];
    const CoordinateType testIntegrationElement =
        testIntegrationElements[testElemPosition];

    // Perform numerical integration
    ResultType sum = 0.;
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const CoordinateType trialWeight =
          trialIntegrationElement * constTrialQuadWeights[trialPoint];
      ResultType partialSum = 0.;
      for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
        const CoordinateType testWeight =
            testIntegrationElement * constTestQuadWeights[testPoint];
        const size_t index = trialPoint * testPointCount * elemPairCount
                             + testPoint * elemPairCount
                             + localElemPairIndex;
        partialSum +=
            kernelValues[index]
            * testBasisValues[testDof * testPointCount + testPoint]
            * trialBasisValues[trialDof * trialPointCount + trialPoint]
            * testWeight;
      }
      sum += partialSum * trialWeight;
    }
    result[i] = sum;
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached(
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

    // Evaluate kernel
    KernelType kernelValues[10 * 10];
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
        Laplace3dSingleLayerPotentialKernel(
            testPointCoo, trialPointCoo,
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
              trialIntegrationElement * constTrialQuadWeights[trialPoint];
          ResultType partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * constTestQuadWeights[testPoint];
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

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached(
    const int elemPairIndexBegin, const unsigned int elemPairCount,
    const unsigned int trialIndexCount,
    const int* __restrict__ testIndices, const int* __restrict__ trialIndices,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const unsigned int testDofCount, BasisFunctionType* testBasisValues,
    const unsigned int trialDofCount, BasisFunctionType* trialBasisValues,
    const unsigned int testElemCount, const unsigned int testVtxCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testVertices,
    const int* testElementCorners,
    const unsigned int trialElemCount, const unsigned int trialVtxCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialVertices,
    const int* trialElementCorners,
    ResultType* __restrict__ result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < elemPairCount) {

    const int coordCount = 3;

    // Determine test and trial element indices
    const int elemPairIndex = elemPairIndexBegin + i;
    const int testElemPosition = testIndices[elemPairIndex / trialIndexCount];
    const int trialElemPosition = trialIndices[elemPairIndex % trialIndexCount];

    // Gather coordinates
    CoordinateType testElemVtx0[coordCount];
    CoordinateType testElemVtx1[coordCount];
    CoordinateType testElemVtx2[coordCount];

    CoordinateType trialElemVtx0[coordCount];
    CoordinateType trialElemVtx1[coordCount];
    CoordinateType trialElemVtx2[coordCount];

    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {

      testElemVtx0[coordIndex] =
          testVertices[testElementCorners[testElemPosition]
                       + coordIndex * testVtxCount];
      testElemVtx1[coordIndex] =
          testVertices[testElementCorners[testElemPosition + testElemCount]
                       + coordIndex * testVtxCount];
      testElemVtx2[coordIndex] =
          testVertices[testElementCorners[testElemPosition + 2 * testElemCount]
                       + coordIndex * testVtxCount];

      trialElemVtx0[coordIndex] =
          trialVertices[trialElementCorners[trialElemPosition]
                        + coordIndex * trialVtxCount];
      trialElemVtx1[coordIndex] =
          trialVertices[trialElementCorners[trialElemPosition + trialElemCount]
                        + coordIndex * trialVtxCount];
      trialElemVtx2[coordIndex] =
          trialVertices[trialElementCorners[trialElemPosition + 2 * trialElemCount]
                        + coordIndex * trialVtxCount];
    }

    // Calculate normals and integration elements
    CoordinateType testElemNormal[coordCount];
    CoordinateType trialElemNormal[coordCount];
    testElemNormal[0] =
        (testElemVtx1[1] - testElemVtx0[1]) * (testElemVtx2[2] - testElemVtx0[2])
      - (testElemVtx1[2] - testElemVtx0[2]) * (testElemVtx2[1] - testElemVtx0[1]);
    testElemNormal[1] =
        (testElemVtx1[2] - testElemVtx0[2]) * (testElemVtx2[0] - testElemVtx0[0])
      - (testElemVtx1[0] - testElemVtx0[0]) * (testElemVtx2[2] - testElemVtx0[2]);
    testElemNormal[2] =
        (testElemVtx1[0] - testElemVtx0[0]) * (testElemVtx2[1] - testElemVtx0[1])
      - (testElemVtx1[1] - testElemVtx0[1]) * (testElemVtx2[0] - testElemVtx0[0]);

    trialElemNormal[0] =
        (trialElemVtx1[1] - trialElemVtx0[1]) * (trialElemVtx2[2] - trialElemVtx0[2])
      - (trialElemVtx1[2] - trialElemVtx0[2]) * (trialElemVtx2[1] - trialElemVtx0[1]);
    trialElemNormal[1] =
        (trialElemVtx1[2] - trialElemVtx0[2]) * (trialElemVtx2[0] - trialElemVtx0[0])
      - (trialElemVtx1[0] - trialElemVtx0[0]) * (trialElemVtx2[2] - trialElemVtx0[2]);
    trialElemNormal[2] =
        (trialElemVtx1[0] - trialElemVtx0[0]) * (trialElemVtx2[1] - trialElemVtx0[1])
      - (trialElemVtx1[1] - trialElemVtx0[1]) * (trialElemVtx2[0] - trialElemVtx0[0]);

    const CoordinateType testIntegrationElement =
        std::sqrt(testElemNormal[0]*testElemNormal[0]
                + testElemNormal[1]*testElemNormal[1]
                + testElemNormal[2]*testElemNormal[2]);
    const CoordinateType trialIntegrationElement =
        std::sqrt(trialElemNormal[0]*trialElemNormal[0]
                + trialElemNormal[1]*trialElemNormal[1]
                + trialElemNormal[2]*trialElemNormal[2]);

    // Calculate global points
    CoordinateType trialElemGlobalPoints[10 * coordCount];
    CoordinateType testElemGlobalPoints[10 * coordCount];
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const CoordinateType ptFun0 = constTrialGeomShapeFun0[trialPoint];
      const CoordinateType ptFun1 = constTrialGeomShapeFun1[trialPoint];
      const CoordinateType ptFun2 = constTrialGeomShapeFun2[trialPoint];
      for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
        trialElemGlobalPoints[coordCount * trialPoint + coordIndex] =
            ptFun0 * trialElemVtx0[coordIndex]
          + ptFun1 * trialElemVtx1[coordIndex]
          + ptFun2 * trialElemVtx2[coordIndex];
      }
    }
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const CoordinateType ptFun0 = constTestGeomShapeFun0[testPoint];
      const CoordinateType ptFun1 = constTestGeomShapeFun1[testPoint];
      const CoordinateType ptFun2 = constTestGeomShapeFun2[testPoint];
      for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
        testElemGlobalPoints[coordCount * testPoint + coordIndex] =
            ptFun0 * testElemVtx0[coordIndex]
          + ptFun1 * testElemVtx1[coordIndex]
          + ptFun2 * testElemVtx2[coordIndex];
      }
    }

    // Evaluate kernel
    KernelType kernelValues[10 * 10];
    CoordinateType trialPointCoo[coordCount], testPointCoo[coordCount];
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
        KernelType kernelValue;
#pragma unroll
        for (size_t coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
          testPointCoo[coordIndex] =
              testElemGlobalPoints[coordCount * testPoint + coordIndex];
          trialPointCoo[coordIndex] =
              trialElemGlobalPoints[coordCount * trialPoint + coordIndex];
        }
        Laplace3dSingleLayerPotentialKernel(
            testPointCoo, trialPointCoo,
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
              trialIntegrationElement * constTrialQuadWeights[trialPoint];
          ResultType partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * constTestQuadWeights[testPoint];
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
