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

#include "cuda.cuh"

#include "../common/common.hpp"
#include "../common/scalar_traits.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct EvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached {

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  QuadData<CoordinateType> testQuadData, trialQuadData;
  BasisFunData<BasisFunctionType> testBasisData, trialBasisData;
  ElemData<CoordinateType> testElemData, trialElemData;
  thrust::device_ptr<ResultType> result;

  EvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData<CoordinateType> _testQuadData,
      const QuadData<CoordinateType> _trialQuadData,
      const BasisFunData<BasisFunctionType> _testBasisData,
      const BasisFunData<BasisFunctionType> _trialBasisData,
      const ElemData<CoordinateType> _testElemData,
      const ElemData<CoordinateType> _trialElemData,
      thrust::device_ptr<ResultType> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testQuadData(_testQuadData), trialQuadData(_trialQuadData),
        testBasisData(_testBasisData), trialBasisData(_trialBasisData),
        testElemData(_testElemData), trialElemData(_trialElemData),
        result (_result) {

    assert(testQuadData.pointCount <= 10);
    assert(trialQuadData.pointCount <= 10);
    assert(testBasisData.dofCount * trialBasisData.dofCount <= 36);
  }

  __host__ __device__
  void operator() (const unsigned int i) {

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    clock_t start_time = clock();

    const CoordinateType testIntegrationElement =
        testElemData.integrationElements[testElemPosition];
    const CoordinateType trialIntegrationElement =
        trialElemData.integrationElements[trialElemPosition];

    const unsigned int activeTestElemCount = testElemData.activeElemCount;
    const unsigned int activeTrialElemCount = trialElemData.activeElemCount;

    const CoordinateType xTestElemNormal = testElemData.normals[testElemPosition];
    const CoordinateType yTestElemNormal = testElemData.normals[testElemPosition+activeTestElemCount];
    const CoordinateType zTestElemNormal = testElemData.normals[testElemPosition+2*activeTestElemCount];

    const CoordinateType xTrialElemNormal = trialElemData.normals[trialElemPosition];
    const CoordinateType yTrialElemNormal = trialElemData.normals[trialElemPosition+activeTrialElemCount];
    const CoordinateType zTrialElemNormal = trialElemData.normals[trialElemPosition+2*activeTrialElemCount];

    clock_t stop_time = clock();

    if (i == 0)
      printf("Gather normals and integration elements: %d clock cycles\n",
          stop_time-start_time);

    start_time = clock();

    CoordinateType xTestGeomData[10], yTestGeomData[10], zTestGeomData[10];
    CoordinateType xTrialGeomData[10], yTrialGeomData[10], zTrialGeomData[10];
    CoordinateType testElemQuadWeights[10], trialElemQuadWeights[10];

    const unsigned int testPointCount = testQuadData.pointCount;
    const unsigned int trialPointCount = trialQuadData.pointCount;

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      xTestGeomData[testPoint] =
          testElemData.geomData[testPoint * activeTestElemCount + testElemPosition];
      yTestGeomData[testPoint] =
          testElemData.geomData[testPoint * activeTestElemCount + testElemPosition
                                + testPointCount * activeTestElemCount];
      zTestGeomData[testPoint] =
          testElemData.geomData[testPoint * activeTestElemCount + testElemPosition
                                + 2 * testPointCount * activeTestElemCount];
      testElemQuadWeights[testPoint] = testQuadData.weights[testPoint];
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      xTrialGeomData[trialPoint] =
          trialElemData.geomData[trialPoint * activeTrialElemCount + trialElemPosition];
      yTrialGeomData[trialPoint] =
          trialElemData.geomData[trialPoint * activeTrialElemCount + trialElemPosition
                                 + trialPointCount * activeTrialElemCount];
      zTrialGeomData[trialPoint] =
          trialElemData.geomData[trialPoint * activeTrialElemCount + trialElemPosition
                                 + 2 * trialPointCount * activeTrialElemCount];
      trialElemQuadWeights[trialPoint] = trialQuadData.weights[trialPoint];
    }

    stop_time = clock();

    if (i == 0)
      printf("Gather global points: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

    ResultType localResult[36];

    const unsigned int testDofCount = testBasisData.dofCount;
    const unsigned int trialDofCount = trialBasisData.dofCount;

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
            const CoordinateType xDist = xTestGeomData[testPoint] - xTrialGeomData[trialPoint];
            const CoordinateType yDist = yTestGeomData[testPoint] - yTrialGeomData[trialPoint];
            const CoordinateType zDist = zTestGeomData[testPoint] - zTrialGeomData[trialPoint];
            const CoordinateType distance = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
            const KernelType kernelValue = static_cast<CoordinateType>(1. / (4. * M_PI)) / sqrt(distance);
            partialSum += kernelValue
//                * testBasisData.basisData[testDof + testDofCount * testPoint]
//                * trialBasisData.basisData[trialDof + trialDofCount * trialPoint]
                * testWeight;
          }
          sum += partialSum * trialWeight;
        }
        localResult[testDof * trialDofCount + trialDof] = sum;
      }
    }

    stop_time = clock();

    if (i == 0)
      printf("Numerical integration: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

    // Copy local result to global device memory
    const unsigned int offset = i * testDofCount * trialDofCount;
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        result[offset + testDof * trialDofCount + trialDof] =
            localResult[testDof * trialDofCount + trialDof];
      }
    }

    stop_time = clock();

    if (i == 0)
      printf("Copy to global memory: %d clock cycles\n", stop_time-start_time);
  }
};

template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct EvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached {

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  QuadData<CoordinateType> testQuadData, trialQuadData;
  BasisFunData<BasisFunctionType> testBasisData, trialBasisData;
  RawGeometryData<CoordinateType> testRawGeometryData, trialRawGeometryData;
  GeomShapeFunData<CoordinateType> testGeomShapeFunData, trialGeomShapeFunData;
  thrust::device_ptr<ResultType> result;

  EvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached(
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
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testQuadData(_testQuadData), trialQuadData(_trialQuadData),
        testBasisData(_testBasisData), trialBasisData(_trialBasisData),
        testRawGeometryData(_testRawGeometryData),
        trialRawGeometryData(_trialRawGeometryData),
        testGeomShapeFunData(_testGeomShapeFunData),
        trialGeomShapeFunData(_trialGeomShapeFunData),
        result (_result) {

    assert(testQuadData.pointCount <= 10);
    assert(trialQuadData.pointCount <= 10);
    assert(testBasisData.dofCount * trialBasisData.dofCount <= 36);
  }

  __device__
  void operator() (const unsigned int i) {

    const int coordCount = 3;

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    clock_t start_time = clock();

    CoordinateType testElemVtx0[coordCount];
    CoordinateType testElemVtx1[coordCount];
    CoordinateType testElemVtx2[coordCount];

    CoordinateType trialElemVtx0[coordCount];
    CoordinateType trialElemVtx1[coordCount];
    CoordinateType trialElemVtx2[coordCount];

    const unsigned int testVtxCount = testRawGeometryData.vtxCount;
    const unsigned int testElemCount = testRawGeometryData.elemCount;

    const unsigned int trialVtxCount = trialRawGeometryData.vtxCount;
    const unsigned int trialElemCount = trialRawGeometryData.elemCount;

    for (int i = 0; i < coordCount; ++i) {

      testElemVtx0[i] = testRawGeometryData.vertices[testRawGeometryData.elementCorners[testElemPosition]+i*testVtxCount];
      testElemVtx1[i] = testRawGeometryData.vertices[testRawGeometryData.elementCorners[testElemPosition+testElemCount]+i*testVtxCount];
      testElemVtx2[i] = testRawGeometryData.vertices[testRawGeometryData.elementCorners[testElemPosition+2*testElemCount]+i*testVtxCount];

      trialElemVtx0[i] = trialRawGeometryData.vertices[trialRawGeometryData.elementCorners[trialElemPosition]+i*trialVtxCount];
      trialElemVtx1[i] = trialRawGeometryData.vertices[trialRawGeometryData.elementCorners[trialElemPosition+trialElemCount]+i*trialVtxCount];
      trialElemVtx2[i] = trialRawGeometryData.vertices[trialRawGeometryData.elementCorners[trialElemPosition+2*trialElemCount]+i*trialVtxCount];
    }

    clock_t stop_time = clock();

    if (i == 0)
      printf("Gather coordinates: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

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

    stop_time = clock();

    if (i == 0)
      printf("Calculate normals and integration elements: %d clock cycles\n",
          stop_time-start_time);

    start_time = clock();

    CoordinateType testGeomData[10 * coordCount], trialGeomData[10 * coordCount];
    CoordinateType testElemQuadWeights[10], trialElemQuadWeights[10];

    const unsigned int testPointCount = testQuadData.pointCount;
    const unsigned int trialPointCount = trialQuadData.pointCount;

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
//      const CoordinateType ptFun0 = testGeomShapeFunData.fun0[testPoint];
//      const CoordinateType ptFun1 = testGeomShapeFunData.fun1[testPoint];
//      const CoordinateType ptFun2 = testGeomShapeFunData.fun2[testPoint];
      const CoordinateType ptFun0 = constTestGeomShapeFun0[testPoint];
      const CoordinateType ptFun1 = constTestGeomShapeFun1[testPoint];
      const CoordinateType ptFun2 = constTestGeomShapeFun2[testPoint];
      for (int i = 0; i < coordCount; ++i) {
        testGeomData[coordCount * testPoint + i] =
            ptFun0 * testElemVtx0[i]
          + ptFun1 * testElemVtx1[i]
          + ptFun2 * testElemVtx2[i];
      }
//      testElemQuadWeights[testPoint] = testQuadData.weights[testPoint];
      testElemQuadWeights[testPoint] = constTestQuadWeights[testPoint];
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
//      const CoordinateType ptFun0 = trialGeomShapeFunData.fun0[trialPoint];
//      const CoordinateType ptFun1 = trialGeomShapeFunData.fun1[trialPoint];
//      const CoordinateType ptFun2 = trialGeomShapeFunData.fun2[trialPoint];
      const CoordinateType ptFun0 = constTrialGeomShapeFun0[trialPoint];
      const CoordinateType ptFun1 = constTrialGeomShapeFun1[trialPoint];
      const CoordinateType ptFun2 = constTrialGeomShapeFun2[trialPoint];
      for (int i = 0; i < coordCount; ++i) {
        trialGeomData[coordCount * trialPoint + i] =
            ptFun0 * trialElemVtx0[i]
          + ptFun1 * trialElemVtx1[i]
          + ptFun2 * trialElemVtx2[i];
      }
//      trialElemQuadWeights[trialPoint] = trialQuadData.weights[trialPoint];
      trialElemQuadWeights[trialPoint] = constTrialQuadWeights[trialPoint];
    }

    stop_time = clock();

    if (i == 0)
      printf("Calculate global points: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

    ResultType localResult[36];

    const unsigned int testDofCount = testBasisData.dofCount;
    const unsigned int trialDofCount = trialBasisData.dofCount;

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

            // evaluate kernel
            CoordinateType sum = 0;
            for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
              CoordinateType diff =
                  testGeomData[coordCount * testPoint + coordIndex]
                - trialGeomData[coordCount * trialPoint + coordIndex];
              sum += diff * diff;
            }
            const KernelType kernelValue = static_cast<CoordinateType>(1. / (4. * M_PI)) / sqrt(sum);

            partialSum += kernelValue
//                * testBasisData.basisData[testDof + testDofCount * testPoint]
//                * trialBasisData.basisData[trialDof + trialDofCount * trialPoint]
                * testWeight;
          }
          sum += partialSum * trialWeight;
        }
        localResult[testDof * trialDofCount + trialDof] = sum;
      }
    }

    stop_time = clock();

    if (i == 0)
      printf("Numerical integration: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

    // Copy local result to global device memory
    const unsigned int offset = i * testDofCount * trialDofCount;
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        result[offset + testDof * trialDofCount + trialDof] =
            localResult[testDof * trialDofCount + trialDof];
      }
    }

    stop_time = clock();

    if (i == 0)
      printf("Copy to global memory: %d clock cycles\n", stop_time-start_time);
  }
};

} // namespace Fiber

#endif
