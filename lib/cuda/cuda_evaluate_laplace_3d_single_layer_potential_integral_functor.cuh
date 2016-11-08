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

#ifndef fiber_cuda_evaluate_laplace_3d_single_layer_potential_integral_functor_cuh
#define fiber_cuda_evaluate_laplace_3d_single_layer_potential_integral_functor_cuh

#include "cuda_evaluate_integral_functor.cuh"
#include "cuda_laplace_3d_single_layer_potential_kernel_functor.hpp"

#include <device_launch_parameters.h>

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
__global__ void
RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached(
    const unsigned int geometryPairCount,
    int* elementPairTestIndexPositions, int* elementPairTrialIndexPositions,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const unsigned int testDofCount, BasisFunctionType* testBasisValues,
    const unsigned int trialDofCount, BasisFunctionType* trialBasisValues,
    const unsigned int testElemCount, const unsigned int testVtxCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* testVertices,
    const int* testElementCorners,
    const unsigned int trialElemCount, const unsigned int trialVtxCount,
    const typename ScalarTraits<BasisFunctionType>::RealType* trialVertices,
    const int* trialElementCorners,
    ResultType* result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  const unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < geometryPairCount) {

//    printf("i=%d\n", i);

    const int coordCount = 3;

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

//    clock_t start_time = clock();

    // Gather coordinates
    CoordinateType testElemVtx0[coordCount];
    CoordinateType testElemVtx1[coordCount];
    CoordinateType testElemVtx2[coordCount];

    CoordinateType trialElemVtx0[coordCount];
    CoordinateType trialElemVtx1[coordCount];
    CoordinateType trialElemVtx2[coordCount];

    for (int coo = 0; coo < coordCount; ++coo) {

      testElemVtx0[coo] =
          testVertices[testElementCorners[testElemPosition]+coo*testVtxCount];
      testElemVtx1[coo] =
          testVertices[testElementCorners[testElemPosition+testElemCount]+coo*testVtxCount];
      testElemVtx2[coo] =
          testVertices[testElementCorners[testElemPosition+2*testElemCount]+coo*testVtxCount];

      trialElemVtx0[coo] =
          trialVertices[trialElementCorners[trialElemPosition]+coo*trialVtxCount];
      trialElemVtx1[coo] =
          trialVertices[trialElementCorners[trialElemPosition+trialElemCount]+coo*trialVtxCount];
      trialElemVtx2[coo] =
          trialVertices[trialElementCorners[trialElemPosition+2*trialElemCount]+coo*trialVtxCount];
    }

//    clock_t stop_time = clock();
//
//    if (i == 0)
//      printf("Gather coordinates: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

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

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Calculate normals and integration elements: %d clock cycles\n",
//          stop_time-start_time);
//
//    start_time = clock();

    // Calculate global points
    CoordinateType testGeomData[10 * coordCount], trialGeomData[10 * coordCount];
    CoordinateType testElemQuadWeights[10], trialElemQuadWeights[10];

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const CoordinateType ptFun0 = constTestGeomShapeFun0[testPoint];
      const CoordinateType ptFun1 = constTestGeomShapeFun1[testPoint];
      const CoordinateType ptFun2 = constTestGeomShapeFun2[testPoint];
      for (int i = 0; i < coordCount; ++i) {
        testGeomData[coordCount * testPoint + i] =
            ptFun0 * testElemVtx0[i]
          + ptFun1 * testElemVtx1[i]
          + ptFun2 * testElemVtx2[i];
      }
      testElemQuadWeights[testPoint] = constTestQuadWeights[testPoint];
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const CoordinateType ptFun0 = constTrialGeomShapeFun0[trialPoint];
      const CoordinateType ptFun1 = constTrialGeomShapeFun1[trialPoint];
      const CoordinateType ptFun2 = constTrialGeomShapeFun2[trialPoint];
      for (int i = 0; i < coordCount; ++i) {
        trialGeomData[coordCount * trialPoint + i] =
            ptFun0 * trialElemVtx0[i]
          + ptFun1 * trialElemVtx1[i]
          + ptFun2 * trialElemVtx2[i];
      }
      trialElemQuadWeights[trialPoint] = constTrialQuadWeights[trialPoint];
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Calculate global points: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

    // Perform numerical integration
    ResultType localResult[36];

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        ResultType sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * trialElemQuadWeights[trialPoint];
          ResultType partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * testElemQuadWeights[testPoint];

            // Evaluate kernel
            KernelType kernelValue;
            CoordinateType testPointCoo[coordCount], trialPointCoo[coordCount];
            for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
              testPointCoo[coordIndex] = testGeomData[coordCount * testPoint + coordIndex];
              trialPointCoo[coordIndex] = trialGeomData[coordCount * trialPoint + coordIndex];
            }
            CudaLaplace3dSingleLayerPotentialKernelFunctor<KernelType>::evaluate(
                testPointCoo, trialPointCoo, kernelValue);

            // Get basis function values
            BasisFunctionType testBasisValue =
                testBasisValues[testDof + testDofCount * testPoint];
            BasisFunctionType trialBasisValue =
                trialBasisValues[trialDof + trialDofCount * trialPoint];

            partialSum += kernelValue * testBasisValue * trialBasisValue * testWeight;
          }
          sum += partialSum * trialWeight;
        }
        localResult[testDof * trialDofCount + trialDof] = sum;
      }
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Numerical integration: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

  // Copy local result to global device memory
  const unsigned int offset = i * testDofCount * trialDofCount;
  for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
    for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
      result[offset + testDof * trialDofCount + trialDof] =
          localResult[testDof * trialDofCount + trialDof];
    }
  }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Copy to global memory: %d clock cycles\n", stop_time-start_time);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached
    : CudaEvaluateIntegralFunctorCached<
      BasisFunctionType, KernelType, ResultType> {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef CudaEvaluateIntegralFunctorCached<
      BasisFunctionType, KernelType, ResultType> Base;

  CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData<CoordinateType> _testQuadData,
      const QuadData<CoordinateType> _trialQuadData,
      const BasisFunData<BasisFunctionType> _testBasisData,
      const BasisFunData<BasisFunctionType> _trialBasisData,
      const ElemData<CoordinateType> _testElemData,
      const ElemData<CoordinateType> _trialElemData,
      thrust::device_ptr<ResultType> _result)
        : Base(
            _elementPairTestIndexPositions, _elementPairTrialIndexPositions,
            _testQuadData, _trialQuadData,
            _testBasisData, _trialBasisData,
            _testElemData, _trialElemData,
            _result) { }

  __device__
  void operator() (const unsigned int i) {

    const int coordCount = 3;

    const int testElemPosition = Base::elementPairTestIndexPositions[i];
    const int trialElemPosition = Base::elementPairTrialIndexPositions[i];

//    clock_t start_time = clock();

    // Gather normals and integration elements
    const CoordinateType testIntegrationElement =
        Base::testElemData.integrationElements[testElemPosition];
    const CoordinateType trialIntegrationElement =
        Base::trialElemData.integrationElements[trialElemPosition];

    const unsigned int activeTestElemCount =
        Base::testElemData.activeElemCount;
    const unsigned int activeTrialElemCount =
        Base::trialElemData.activeElemCount;

//    CoordinateType testElemNormal[coordCount];
//    CoordinateType trialElemNormal[coordCount];
//
//    testElemNormal[0] = Base::testElemData.normals[testElemPosition];
//    testElemNormal[1] = Base::testElemData.normals[testElemPosition+activeTestElemCount];
//    testElemNormal[2] = Base::testElemData.normals[testElemPosition+2*activeTestElemCount];
//
//    trialElemNormal[0] = Base::trialElemData.normals[trialElemPosition];
//    trialElemNormal[1] = Base::trialElemData.normals[trialElemPosition+activeTrialElemCount];
//    trialElemNormal[2] = Base::trialElemData.normals[trialElemPosition+2*activeTrialElemCount];

//    clock_t stop_time = clock();
//
//    if (i == 0)
//      printf("Gather normals and integration elements: %d clock cycles\n",
//          stop_time-start_time);
//
//    start_time = clock();

    // Gather globals points
    CoordinateType testGeomData[10 * coordCount], trialGeomData[10 * coordCount];
    CoordinateType testElemQuadWeights[10], trialElemQuadWeights[10];

    const unsigned int testPointCount = Base::testQuadData.pointCount;
    const unsigned int trialPointCount = Base::trialQuadData.pointCount;

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      for (int coo = 0; coo < coordCount; ++coo) {
        testGeomData[coordCount * testPoint + coo] =
            Base::testElemData.geomData[
                testPoint * activeTestElemCount
                + testElemPosition
                + coo * testPointCount * activeTestElemCount];
      }
      testElemQuadWeights[testPoint] = constTestQuadWeights[testPoint];
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      for (int coo = 0; coo < coordCount; ++coo) {
        trialGeomData[coordCount * trialPoint + coo] =
            Base::trialElemData.geomData[
                trialPoint * activeTrialElemCount
                + trialElemPosition
                + coo * trialPointCount * activeTrialElemCount];
      }
      trialElemQuadWeights[trialPoint] = constTrialQuadWeights[trialPoint];
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Gather global points: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

    // Perform numerical integration
    ResultType localResult[36];

    const unsigned int testDofCount = Base::testBasisData.dofCount;
    const unsigned int trialDofCount = Base::trialBasisData.dofCount;

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        ResultType sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * trialElemQuadWeights[trialPoint];
          ResultType partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * testElemQuadWeights[testPoint];

            // Evaluate kernel
            KernelType kernelValue;
            CoordinateType testPointCoo[coordCount], trialPointCoo[coordCount];
            for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
              testPointCoo[coordIndex] = testGeomData[coordCount * testPoint + coordIndex];
              trialPointCoo[coordIndex] = trialGeomData[coordCount * trialPoint + coordIndex];
            }
            CudaLaplace3dSingleLayerPotentialKernelFunctor<KernelType>::evaluate(
                testPointCoo, trialPointCoo, kernelValue);

            // Get basis function values
            BasisFunctionType testBasisValue =
                Base::testBasisData.values[testDof + testDofCount * testPoint];
            BasisFunctionType trialBasisValue =
                Base::trialBasisData.values[trialDof + trialDofCount * trialPoint];

            partialSum += kernelValue * testBasisValue * trialBasisValue * testWeight;
          }
          sum += partialSum * trialWeight;
        }
        localResult[testDof * trialDofCount + trialDof] = sum;
      }
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Numerical integration: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

    // Copy local result to global device memory
    const unsigned int offset = i * testDofCount * trialDofCount;
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        Base::result[offset + testDof * trialDofCount + trialDof] =
            localResult[testDof * trialDofCount + trialDof];
      }
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Copy to global memory: %d clock cycles\n", stop_time-start_time);
  }
};

template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached
    : CudaEvaluateIntegralFunctorNonCached<
      BasisFunctionType, KernelType, ResultType> {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef CudaEvaluateIntegralFunctorNonCached<
      BasisFunctionType, KernelType, ResultType> Base;

  CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData<CoordinateType> _testQuadData,
      const QuadData<CoordinateType> _trialQuadData,
      const BasisFunData<BasisFunctionType> _testBasisData,
      const BasisFunData<BasisFunctionType> _trialBasisData,
      const RawGeometryData<CoordinateType> _testRawGeometryData,
      const RawGeometryData<CoordinateType> _trialRawGeometryData,
      const GeomShapeFunData<CoordinateType> _testGeomShapeFunData,
      const GeomShapeFunData<CoordinateType> _trialGeomShapeFunData,
      thrust::device_ptr<ResultType> _result)
      : Base(
          _elementPairTestIndexPositions, _elementPairTrialIndexPositions,
          _testQuadData, _trialQuadData,
          _testBasisData, _trialBasisData,
          _testRawGeometryData, _trialRawGeometryData,
          _testGeomShapeFunData, _trialGeomShapeFunData,
          _result) { }

  __device__
  void operator() (const unsigned int i) {

    const int coordCount = 3;

    const int testElemPosition = Base::elementPairTestIndexPositions[i];
    const int trialElemPosition = Base::elementPairTrialIndexPositions[i];

//    clock_t start_time = clock();

    // Gather coordinates
    CoordinateType testElemVtx0[coordCount];
    CoordinateType testElemVtx1[coordCount];
    CoordinateType testElemVtx2[coordCount];

    CoordinateType trialElemVtx0[coordCount];
    CoordinateType trialElemVtx1[coordCount];
    CoordinateType trialElemVtx2[coordCount];

    const unsigned int testVtxCount = Base::testRawGeometryData.vtxCount;
    const unsigned int testElemCount = Base::testRawGeometryData.elemCount;

    const unsigned int trialVtxCount = Base::trialRawGeometryData.vtxCount;
    const unsigned int trialElemCount = Base::trialRawGeometryData.elemCount;

    for (int coo = 0; coo < coordCount; ++coo) {

      testElemVtx0[coo] =
          Base::testRawGeometryData.vertices[Base::testRawGeometryData.elementCorners[testElemPosition]+coo*testVtxCount];
      testElemVtx1[coo] =
          Base::testRawGeometryData.vertices[Base::testRawGeometryData.elementCorners[testElemPosition+testElemCount]+coo*testVtxCount];
      testElemVtx2[coo] =
          Base::testRawGeometryData.vertices[Base::testRawGeometryData.elementCorners[testElemPosition+2*testElemCount]+coo*testVtxCount];

      trialElemVtx0[coo] =
          Base::trialRawGeometryData.vertices[Base::trialRawGeometryData.elementCorners[trialElemPosition]+coo*trialVtxCount];
      trialElemVtx1[coo] =
          Base::trialRawGeometryData.vertices[Base::trialRawGeometryData.elementCorners[trialElemPosition+trialElemCount]+coo*trialVtxCount];
      trialElemVtx2[coo] =
          Base::trialRawGeometryData.vertices[Base::trialRawGeometryData.elementCorners[trialElemPosition+2*trialElemCount]+coo*trialVtxCount];
    }

//    clock_t stop_time = clock();
//
//    if (i == 0)
//      printf("Gather coordinates: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

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

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Calculate normals and integration elements: %d clock cycles\n",
//          stop_time-start_time);
//
//    start_time = clock();

    // Calculate global points
    CoordinateType testGeomData[10 * coordCount], trialGeomData[10 * coordCount];
    CoordinateType testElemQuadWeights[10], trialElemQuadWeights[10];

    const unsigned int testPointCount = Base::testQuadData.pointCount;
    const unsigned int trialPointCount = Base::trialQuadData.pointCount;

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const CoordinateType ptFun0 = constTestGeomShapeFun0[testPoint];
      const CoordinateType ptFun1 = constTestGeomShapeFun1[testPoint];
      const CoordinateType ptFun2 = constTestGeomShapeFun2[testPoint];
      for (int coo = 0; coo < coordCount; ++coo) {
        testGeomData[coordCount * testPoint + coo] =
            ptFun0 * testElemVtx0[coo]
          + ptFun1 * testElemVtx1[coo]
          + ptFun2 * testElemVtx2[coo];
      }
      testElemQuadWeights[testPoint] = constTestQuadWeights[testPoint];
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const CoordinateType ptFun0 = constTrialGeomShapeFun0[trialPoint];
      const CoordinateType ptFun1 = constTrialGeomShapeFun1[trialPoint];
      const CoordinateType ptFun2 = constTrialGeomShapeFun2[trialPoint];
      for (int coo = 0; coo < coordCount; ++coo) {
        trialGeomData[coordCount * trialPoint + coo] =
            ptFun0 * trialElemVtx0[coo]
          + ptFun1 * trialElemVtx1[coo]
          + ptFun2 * trialElemVtx2[coo];
      }
      trialElemQuadWeights[trialPoint] = constTrialQuadWeights[trialPoint];
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Calculate global points: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

    // Perform numerical integration
    ResultType localResult[36];

    const unsigned int testDofCount = Base::testBasisData.dofCount;
    const unsigned int trialDofCount = Base::trialBasisData.dofCount;

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        ResultType sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const CoordinateType trialWeight =
              trialIntegrationElement * trialElemQuadWeights[trialPoint];
          ResultType partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const CoordinateType testWeight =
                testIntegrationElement * testElemQuadWeights[testPoint];

            // Evaluate kernel
            KernelType kernelValue;
            CoordinateType testPointCoo[coordCount], trialPointCoo[coordCount];
            for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
              testPointCoo[coordIndex] = testGeomData[coordCount * testPoint + coordIndex];
              trialPointCoo[coordIndex] = trialGeomData[coordCount * trialPoint + coordIndex];
            }
            CudaLaplace3dSingleLayerPotentialKernelFunctor<KernelType>::evaluate(
                testPointCoo, trialPointCoo, kernelValue);

            // Get basis function values
            BasisFunctionType testBasisValue =
                Base::testBasisData.values[testDof + testDofCount * testPoint];
            BasisFunctionType trialBasisValue =
                Base::trialBasisData.values[trialDof + trialDofCount * trialPoint];

            partialSum += kernelValue * testBasisValue * trialBasisValue * testWeight;
          }
          sum += partialSum * trialWeight;
        }
        localResult[testDof * trialDofCount + trialDof] = sum;
      }
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Numerical integration: %d clock cycles\n", stop_time-start_time);
//
//    start_time = clock();

    // Copy local result to global device memory
    const unsigned int offset = i * testDofCount * trialDofCount;
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        Base::result[offset + testDof * trialDofCount + trialDof] =
            localResult[testDof * trialDofCount + trialDof];
      }
    }

//    stop_time = clock();
//
//    if (i == 0)
//      printf("Copy to global memory: %d clock cycles\n", stop_time-start_time);
  }
};

} // namespace Fiber

#endif
