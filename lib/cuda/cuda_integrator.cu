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

#include "cuda_integrator.hpp"
#include "cuda_grid.hpp"

#include "../common/not_implemented_error.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/shapeset.hpp"
#include "../fiber/basis_data.hpp"
#include "../fiber/scalar_traits.hpp"

#include <complex>
#include <chrono>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

namespace Fiber {

template <typename ValueType>
__host__ __device__ void evaluateLaplace3dSingleLayerPotentialKernel(
    double* testGlobalData, double* trialGlobalData, double* result) {

  const int coordCount = 3;

  ValueType sum = 0;
  for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
    ValueType diff =
        testGlobalData[coordIndex] - trialGlobalData[coordIndex];
    sum += diff * diff;
  }
  *result = static_cast<double>(1. / (4. * M_PI)) / sqrt(sum);
}

//template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct evaluateIntegralFunctor {

  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  unsigned int testPointCount;
  unsigned int trialPointCount;
  unsigned int testDofCount;
  unsigned int trialDofCount;
  thrust::device_ptr<double> testQuadWeights;
  thrust::device_ptr<double> trialQuadWeights;
  thrust::device_ptr<double> testBasisData;
  thrust::device_ptr<double> trialBasisData;
  thrust::device_ptr<double> testGeomData;
  thrust::device_ptr<double> trialGeomData;
  unsigned int activeTestElemCount;
  unsigned int activeTrialElemCount;
  thrust::device_ptr<const double> testNormals;
  thrust::device_ptr<const double> trialNormals;
  thrust::device_ptr<const double> testIntegrationElements;
  thrust::device_ptr<const double> trialIntegrationElements;
  thrust::device_ptr<double> result;

  evaluateIntegralFunctor(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const unsigned int _testPointCount,
      const unsigned int _trialPointCount,
      const unsigned int _testDofCount,
      const unsigned int _trialDofCount,
      thrust::device_ptr<double> _testQuadWeights,
      thrust::device_ptr<double> _trialQuadWeights,
      thrust::device_ptr<double> _testBasisData,
      thrust::device_ptr<double> _trialBasisData,
      thrust::device_ptr<double> _testGeomData,
      thrust::device_ptr<double> _trialGeomData,
      const unsigned int _activeTestElemCount,
      const unsigned int _activeTrialElemCount,
      thrust::device_ptr<const double> _testNormals,
      thrust::device_ptr<const double> _trialNormals,
      thrust::device_ptr<const double> _testIntegrationElements,
      thrust::device_ptr<const double> _trialIntegrationElements,
      thrust::device_ptr<double> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testPointCount(_testPointCount),
        trialPointCount(_trialPointCount),
        testDofCount(_testDofCount),
        trialDofCount(_trialDofCount),
        testQuadWeights(_testQuadWeights),
        trialQuadWeights(_trialQuadWeights),
        testBasisData(_testBasisData),
        trialBasisData(_trialBasisData),
        testGeomData(_testGeomData),
        trialGeomData(_trialGeomData),
        activeTestElemCount(_activeTestElemCount),
        activeTrialElemCount(_activeTrialElemCount),
        testNormals(_testNormals),
        trialNormals(_trialNormals),
        testIntegrationElements(_testIntegrationElements),
        trialIntegrationElements(_trialIntegrationElements),
        result (_result) {

    assert(testPointCount <= 10);
    assert(trialPointCount <= 10);
    assert(testDofCount * trialDofCount <= 36);
  }

  __host__ __device__
  void operator() (const unsigned int i) {

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    clock_t start_time = clock();

    const double testIntegrationElement = testIntegrationElements[testElemPosition];
    const double trialIntegrationElement = trialIntegrationElements[trialElemPosition];

    const double xTestElemNormal = testNormals[testElemPosition];
    const double yTestElemNormal = testNormals[testElemPosition+activeTestElemCount];
    const double zTestElemNormal = testNormals[testElemPosition+2*activeTestElemCount];

    const double xTrialElemNormal = trialNormals[trialElemPosition];
    const double yTrialElemNormal = trialNormals[trialElemPosition+activeTrialElemCount];
    const double zTrialElemNormal = trialNormals[trialElemPosition+2*activeTrialElemCount];

    clock_t stop_time = clock();

    if (i == 0)
      printf("Gather normals and integration elements: %d clock cycles\n",
          stop_time-start_time);

    start_time = clock();

//    double* xTestGeomData = new double[testPointCount];
//    double* yTestGeomData = new double[testPointCount];
//    double* zTestGeomData = new double[testPointCount];
//    double* xTrialGeomData = new double[trialPointCount];
//    double* yTrialGeomData = new double[trialPointCount];
//    double* zTrialGeomData = new double[trialPointCount];

    double xTestGeomData[10];
    double yTestGeomData[10];
    double zTestGeomData[10];
    double xTrialGeomData[10];
    double yTrialGeomData[10];
    double zTrialGeomData[10];
    double testElemQuadWeights[10];
    double trialElemQuadWeights[10];

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      xTestGeomData[testPoint] =
          testGeomData[testPoint * activeTestElemCount + testElemPosition];
      yTestGeomData[testPoint] =
          testGeomData[testPoint * activeTestElemCount + testElemPosition
                       + testPointCount * activeTestElemCount];
      zTestGeomData[testPoint] =
          testGeomData[testPoint * activeTestElemCount + testElemPosition
                       + 2 * testPointCount * activeTestElemCount];
      testElemQuadWeights[testPoint] = testQuadWeights[testPoint];
    }

    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      xTrialGeomData[trialPoint] =
          trialGeomData[trialPoint * activeTrialElemCount + trialElemPosition];
      yTrialGeomData[trialPoint] =
          trialGeomData[trialPoint * activeTrialElemCount + trialElemPosition
                        + trialPointCount * activeTrialElemCount];
      zTrialGeomData[trialPoint] =
          trialGeomData[trialPoint * activeTrialElemCount + trialElemPosition
                        + 2 * trialPointCount * activeTrialElemCount];
      trialElemQuadWeights[trialPoint] = trialQuadWeights[trialPoint];
    }

    stop_time = clock();

    if (i == 0)
      printf("Gather global points: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

//    double* localResult = new double[testDofCount * trialDofCount];
    double localResult[36];
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        double sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const double trialWeight =
              trialIntegrationElement * trialElemQuadWeights[trialPoint];
          double partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const double testWeight =
                testIntegrationElement * testElemQuadWeights[testPoint];
            // TODO: Check impact of function call on performance
            const double xDist = xTestGeomData[testPoint] - xTrialGeomData[trialPoint];
            const double yDist = yTestGeomData[testPoint] - yTrialGeomData[trialPoint];
            const double zDist = zTestGeomData[testPoint] - zTrialGeomData[trialPoint];
            const double distance = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
            const double kernelValue = 1.0 / (4.0 * M_PI * distance);
//            double* kernelValue;
//            evaluateLaplace3dSingleLayerPotentialKernel(
//                xTestGeomData, TrialGeomData, kernelValue);
            partialSum += kernelValue *
//                          testBasisData[testDof + testDofCount * testPoint] *
//                          trialBasisData[trialDof + trialDofCount * trialPoint] *
                          testWeight;
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

    // Free local arrays
//    delete[] xTestGeomData;
//    delete[] yTestGeomData;
//    delete[] zTestGeomData;
//    delete[] xTrialGeomData;
//    delete[] yTrialGeomData;
//    delete[] zTrialGeomData;
//    delete[] localResult;
  }
};

//template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct evaluateIntegralFunctorNonCached {

  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  unsigned int testPointCount;
  unsigned int trialPointCount;
  unsigned int testDofCount;
  unsigned int trialDofCount;
  thrust::device_ptr<double> testQuadWeights;
  thrust::device_ptr<double> trialQuadWeights;
  thrust::device_ptr<double> testBasisData;
  thrust::device_ptr<double> trialBasisData;
  unsigned int activeTestElemCount;
  unsigned int activeTrialElemCount;
  thrust::device_ptr<const double> testVtx0x;
  thrust::device_ptr<const double> testVtx0y;
  thrust::device_ptr<const double> testVtx0z;
  thrust::device_ptr<const double> testVtx1x;
  thrust::device_ptr<const double> testVtx1y;
  thrust::device_ptr<const double> testVtx1z;
  thrust::device_ptr<const double> testVtx2x;
  thrust::device_ptr<const double> testVtx2y;
  thrust::device_ptr<const double> testVtx2z;
  thrust::device_ptr<const double> trialVtx0x;
  thrust::device_ptr<const double> trialVtx0y;
  thrust::device_ptr<const double> trialVtx0z;
  thrust::device_ptr<const double> trialVtx1x;
  thrust::device_ptr<const double> trialVtx1y;
  thrust::device_ptr<const double> trialVtx1z;
  thrust::device_ptr<const double> trialVtx2x;
  thrust::device_ptr<const double> trialVtx2y;
  thrust::device_ptr<const double> trialVtx2z;
  thrust::device_ptr<const double> testFun0;
  thrust::device_ptr<const double> testFun1;
  thrust::device_ptr<const double> testFun2;
  thrust::device_ptr<const double> trialFun0;
  thrust::device_ptr<const double> trialFun1;
  thrust::device_ptr<const double> trialFun2;
  thrust::device_ptr<double> result;

  evaluateIntegralFunctorNonCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const unsigned int _testPointCount,
      const unsigned int _trialPointCount,
      const unsigned int _testDofCount,
      const unsigned int _trialDofCount,
      thrust::device_ptr<double> _testQuadWeights,
      thrust::device_ptr<double> _trialQuadWeights,
      thrust::device_ptr<double> _testBasisData,
      thrust::device_ptr<double> _trialBasisData,
      const unsigned int _activeTestElemCount,
      const unsigned int _activeTrialElemCount,
      thrust::device_ptr<const double> _testVtx0x,
      thrust::device_ptr<const double> _testVtx0y,
      thrust::device_ptr<const double> _testVtx0z,
      thrust::device_ptr<const double> _testVtx1x,
      thrust::device_ptr<const double> _testVtx1y,
      thrust::device_ptr<const double> _testVtx1z,
      thrust::device_ptr<const double> _testVtx2x,
      thrust::device_ptr<const double> _testVtx2y,
      thrust::device_ptr<const double> _testVtx2z,
      thrust::device_ptr<const double> _trialVtx0x,
      thrust::device_ptr<const double> _trialVtx0y,
      thrust::device_ptr<const double> _trialVtx0z,
      thrust::device_ptr<const double> _trialVtx1x,
      thrust::device_ptr<const double> _trialVtx1y,
      thrust::device_ptr<const double> _trialVtx1z,
      thrust::device_ptr<const double> _trialVtx2x,
      thrust::device_ptr<const double> _trialVtx2y,
      thrust::device_ptr<const double> _trialVtx2z,
      thrust::device_ptr<const double> _testFun0,
      thrust::device_ptr<const double> _testFun1,
      thrust::device_ptr<const double> _testFun2,
      thrust::device_ptr<const double> _trialFun0,
      thrust::device_ptr<const double> _trialFun1,
      thrust::device_ptr<const double> _trialFun2,
      thrust::device_ptr<double> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testPointCount(_testPointCount),
        trialPointCount(_trialPointCount),
        testDofCount(_testDofCount),
        trialDofCount(_trialDofCount),
        testQuadWeights(_testQuadWeights),
        trialQuadWeights(_trialQuadWeights),
        testBasisData(_testBasisData),
        trialBasisData(_trialBasisData),
        activeTestElemCount(_activeTestElemCount),
        activeTrialElemCount(_activeTrialElemCount),
        testVtx0x(_testVtx0x),
        testVtx0y(_testVtx0y),
        testVtx0z(_testVtx0z),
        testVtx1x(_testVtx1x),
        testVtx1y(_testVtx1y),
        testVtx1z(_testVtx1z),
        testVtx2x(_testVtx2x),
        testVtx2y(_testVtx2y),
        testVtx2z(_testVtx2z),
        trialVtx0x(_trialVtx0x),
        trialVtx0y(_trialVtx0y),
        trialVtx0z(_trialVtx0z),
        trialVtx1x(_trialVtx1x),
        trialVtx1y(_trialVtx1y),
        trialVtx1z(_trialVtx1z),
        trialVtx2x(_trialVtx2x),
        trialVtx2y(_trialVtx2y),
        trialVtx2z(_trialVtx2z),
        testFun0(_testFun0),
        testFun1(_testFun1),
        testFun2(_testFun2),
        trialFun0(_trialFun0),
        trialFun1(_trialFun1),
        trialFun2(_trialFun2),
        result (_result) {

    assert(testPointCount <= 10);
    assert(trialPointCount <= 10);
    assert(testDofCount * trialDofCount <= 36);
  }

  __host__ __device__
  void operator() (const unsigned int i) {

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    clock_t start_time = clock();

    const double testElemVtx0x = testVtx0x[testElemPosition];
    const double testElemVtx0y = testVtx0y[testElemPosition];
    const double testElemVtx0z = testVtx0z[testElemPosition];
    const double testElemVtx1x = testVtx1x[testElemPosition];
    const double testElemVtx1y = testVtx1y[testElemPosition];
    const double testElemVtx1z = testVtx1z[testElemPosition];
    const double testElemVtx2x = testVtx2x[testElemPosition];
    const double testElemVtx2y = testVtx2y[testElemPosition];
    const double testElemVtx2z = testVtx2z[testElemPosition];

    const double trialElemVtx0x = trialVtx0x[trialElemPosition];
    const double trialElemVtx0y = trialVtx0y[trialElemPosition];
    const double trialElemVtx0z = trialVtx0z[trialElemPosition];
    const double trialElemVtx1x = trialVtx1x[trialElemPosition];
    const double trialElemVtx1y = trialVtx1y[trialElemPosition];
    const double trialElemVtx1z = trialVtx1z[trialElemPosition];
    const double trialElemVtx2x = trialVtx2x[trialElemPosition];
    const double trialElemVtx2y = trialVtx2y[trialElemPosition];
    const double trialElemVtx2z = trialVtx2z[trialElemPosition];

    clock_t stop_time = clock();

    if (i == 0)
      printf("Gather coordinates: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

    double xTestElemNormal =
        (testElemVtx1y - testElemVtx0y) * (testElemVtx2z - testElemVtx0z)
      - (testElemVtx1z - testElemVtx0z) * (testElemVtx2y - testElemVtx0y);
    double yTestElemNormal =
        (testElemVtx1z - testElemVtx0z) * (testElemVtx2x - testElemVtx0x)
      - (testElemVtx1x - testElemVtx0x) * (testElemVtx2z - testElemVtx0z);
    double zTestElemNormal =
        (testElemVtx1x - testElemVtx0x) * (testElemVtx2y - testElemVtx0y)
      - (testElemVtx1y - testElemVtx0y) * (testElemVtx2x - testElemVtx0x);

    double xTrialElemNormal =
        (trialElemVtx1y - trialElemVtx0y) * (trialElemVtx2z - trialElemVtx0z)
      - (trialElemVtx1z - trialElemVtx0z) * (trialElemVtx2y - trialElemVtx0y);
    double yTrialElemNormal =
        (trialElemVtx1z - trialElemVtx0z) * (trialElemVtx2x - trialElemVtx0x)
      - (trialElemVtx1x - trialElemVtx0x) * (trialElemVtx2z - trialElemVtx0z);
    double zTrialElemNormal =
        (trialElemVtx1x - trialElemVtx0x) * (trialElemVtx2y - trialElemVtx0y)
      - (trialElemVtx1y - trialElemVtx0y) * (trialElemVtx2x - trialElemVtx0x);

    const double testIntegrationElement =
        std::sqrt(xTestElemNormal*xTestElemNormal
                + yTestElemNormal*yTestElemNormal
                + zTestElemNormal*zTestElemNormal);
    const double trialIntegrationElement =
        std::sqrt(xTrialElemNormal*xTrialElemNormal
                + yTrialElemNormal*yTrialElemNormal
                + zTrialElemNormal*zTrialElemNormal);

    xTestElemNormal /= testIntegrationElement;
    yTestElemNormal /= testIntegrationElement;
    zTestElemNormal /= testIntegrationElement;
    xTrialElemNormal /= trialIntegrationElement;
    yTrialElemNormal /= trialIntegrationElement;
    zTrialElemNormal /= trialIntegrationElement;

    stop_time = clock();

    if (i == 0)
      printf("Calculate normals and integration elements: %d clock cycles\n",
          stop_time-start_time);

    start_time = clock();

//    double* xTestGeomData = new double[testPointCount];
//    double* yTestGeomData = new double[testPointCount];
//    double* zTestGeomData = new double[testPointCount];
//    double* xTrialGeomData = new double[trialPointCount];
//    double* yTrialGeomData = new double[trialPointCount];
//    double* zTrialGeomData = new double[trialPointCount];

    double xTestGeomData[10];
    double yTestGeomData[10];
    double zTestGeomData[10];
    double xTrialGeomData[10];
    double yTrialGeomData[10];
    double zTrialGeomData[10];
    double testElemQuadWeights[10];
    double trialElemQuadWeights[10];

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {

      const double ptFun0 = testFun0[testPoint];
      const double ptFun1 = testFun1[testPoint];
      const double ptFun2 = testFun2[testPoint];

      xTestGeomData[testPoint] = ptFun0 * testElemVtx0x
                               + ptFun1 * testElemVtx1x
                               + ptFun2 * testElemVtx2x;;
      yTestGeomData[testPoint] = ptFun0 * testElemVtx0y
                               + ptFun1 * testElemVtx1y
                               + ptFun2 * testElemVtx2y;;
      zTestGeomData[testPoint] = ptFun0 * testElemVtx0z
                               + ptFun1 * testElemVtx1z
                               + ptFun2 * testElemVtx2z;

      testElemQuadWeights[testPoint] = testQuadWeights[testPoint];
    }

    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {

      const double ptFun0 = trialFun0[trialPoint];
      const double ptFun1 = trialFun1[trialPoint];
      const double ptFun2 = trialFun2[trialPoint];

      xTrialGeomData[trialPoint] = ptFun0 * trialElemVtx0x
                                 + ptFun1 * trialElemVtx1x
                                 + ptFun2 * trialElemVtx2x;;
      yTrialGeomData[trialPoint] = ptFun0 * trialElemVtx0y
                                 + ptFun1 * trialElemVtx1y
                                 + ptFun2 * trialElemVtx2y;;
      zTrialGeomData[trialPoint] = ptFun0 * trialElemVtx0z
                                 + ptFun1 * trialElemVtx1z
                                 + ptFun2 * trialElemVtx2z;

      trialElemQuadWeights[trialPoint] = trialQuadWeights[trialPoint];
    }

    stop_time = clock();

    if (i == 0)
      printf("Calculate global points: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

//    double* localResult = new double[testDofCount * trialDofCount];
    double localResult[36];
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        double sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const double trialWeight =
              trialIntegrationElement * trialElemQuadWeights[trialPoint];
          double partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const double testWeight =
                testIntegrationElement * testElemQuadWeights[testPoint];
            // TODO: Check impact of function call on performance
            const double xDist = xTestGeomData[testPoint] - xTrialGeomData[trialPoint];
            const double yDist = yTestGeomData[testPoint] - yTrialGeomData[trialPoint];
            const double zDist = zTestGeomData[testPoint] - zTrialGeomData[trialPoint];
            const double distance = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
            const double kernelValue = 1.0 / (4.0 * M_PI * distance);
//            double* kernelValue;
//            evaluateLaplace3dSingleLayerPotentialKernel(
//                xTestGeomData, TrialGeomData, kernelValue);
            partialSum += kernelValue *
//                          testBasisData[testDof + testDofCount * testPoint] *
//                          trialBasisData[trialDof + trialDofCount * trialPoint] *
                          testWeight;
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

    // Free local arrays
//    delete[] xTestGeomData;
//    delete[] yTestGeomData;
//    delete[] zTestGeomData;
//    delete[] xTrialGeomData;
//    delete[] yTrialGeomData;
//    delete[] zTrialGeomData;
//    delete[] localResult;
  }
};

struct geomShapeFunFunctor {

  __host__ __device__
  thrust::tuple<double, double, double> operator()(
    const thrust::tuple<double, double>& localPointCoo) const {

    const double r = thrust::get<0>(localPointCoo);
    const double s = thrust::get<1>(localPointCoo);

    const double fun0 = 1.0 - r - s;
    const double fun1 = r;
    const double fun2 = s;

    return thrust::make_tuple(fun0, fun1, fun2);
  }
};

template <typename BasisFunctionType, typename ResultType>
CudaIntegrator<BasisFunctionType, ResultType>::CudaIntegrator(
    const Matrix<double> &localTestQuadPoints,
    const Matrix<double> &localTrialQuadPoints,
    const std::vector<double> &testQuadWeights,
    const std::vector<double> &trialQuadWeights,
    shared_ptr<Bempp::CudaGrid> testGrid,
    shared_ptr<Bempp::CudaGrid> trialGrid,
    bool cacheElemData)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testQuadWeights(testQuadWeights), m_trialQuadWeights(trialQuadWeights),
      m_testGrid(testGrid), m_trialGrid(trialGrid),
      m_cacheElemData(cacheElemData) {

  if (localTestQuadPoints.cols() != testQuadWeights.size())
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "numbers of test points and weights do not match");
  if (localTrialQuadPoints.cols() != trialQuadWeights.size())
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "numbers of trial points and weights do not match");
}

template <typename BasisFunctionType, typename ResultType>
CudaIntegrator<BasisFunctionType, ResultType>::~CudaIntegrator() { }

template <typename BasisFunctionType, typename ResultType>
void CudaIntegrator<BasisFunctionType, ResultType>::integrate(
    const std::vector<int> &elementPairTestIndices,
    const std::vector<int> &elementPairTrialIndices,
    const std::vector<int> &testDeviceElemIndices,
    const std::vector<int> &trialDeviceElemIndices,
    const Shapeset<BasisFunctionType> &testShapeset,
    const Shapeset<BasisFunctionType> &trialShapeset,
    std::vector<Matrix<ResultType>*> &result) {

  std::cout << "Hello, this is CudaIntegrator::integrate()!" << std::endl;

  const int testPointCount = m_localTestQuadPoints.cols();
  const int trialPointCount = m_localTrialQuadPoints.cols();
  const int geometryPairCount = elementPairTestIndices.size();

  if (elementPairTestIndices.size() != elementPairTrialIndices.size())
    throw std::invalid_argument(
        "CudaIntegrator::integrate(): "
        "arrays 'elementPairTestIndices' and 'elementPairTrialIndices' must "
        "have the same number of elements");

  if (result.size() != geometryPairCount)
    throw std::invalid_argument(
        "CudaIntegrator::integrate(): "
        "arrays 'result' and 'elementPairIndices' must have the same number "
        "of elements");

  if (testPointCount == 0 || trialPointCount == 0 || geometryPairCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryPairCount != 0, set elements of result to 0.

  const int testDofCount = testShapeset.size();
  const int trialDofCount = trialShapeset.size();

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->resize(testDofCount, trialDofCount);
  }

  // Evaluate shapesets
  BasisData<double> testBasisData, trialBasisData;
  size_t testBasisDeps = 0, trialBasisDeps = 0;
//  testShapeset.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
//                        testBasisData);
//  trialShapeset.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
//                         trialBasisData);
  thrust::device_vector<double> d_testBasisData(
      testBasisData.values.begin(), testBasisData.values.end());
  thrust::device_vector<double> d_trialBasisData(
      trialBasisData.values.begin(), trialBasisData.values.end());

  // Copy numerical quadrature weights to device memory
  thrust::device_vector<double> d_testQuadWeights(m_testQuadWeights);
  thrust::device_vector<double> d_trialQuadWeights(m_trialQuadWeights);

  // Copy element pair indices to device memory
//  thrust::device_ptr<int> d_elementPairTestIndicesPtr;
//  thrust::device_ptr<int> d_elementPairTrialIndicesPtr;
//  d_elementPairTestIndicesPtr = thrust::device_malloc(geometryPairCount);
//  thrust::copy(elementPairTestIndices.begin(), elementPairTestIndices.end(),
//      d_elementPairTestIndicesPtr);
  thrust::device_vector<int> d_elementPairTestIndices(elementPairTestIndices);
  thrust::device_vector<int> d_elementPairTrialIndices(elementPairTrialIndices);

  // Replace element pair indices by their positions in activeElemIndices
  // vectors (not necessary in case of ALL_ELEMS active)
  if (testDeviceElemIndices[0] != -1) {
    for (int testElemPosition = 0; testElemPosition < testDeviceElemIndices.size(); ++testElemPosition) {
      const int activeTestElemIndex = testDeviceElemIndices[testElemPosition];
      thrust::replace(d_elementPairTestIndices.begin(), d_elementPairTestIndices.end(),
          activeTestElemIndex, testElemPosition);
    }
  }
  if (trialDeviceElemIndices[0] != -1) {
    for (int trialElemPosition = 0; trialElemPosition < trialDeviceElemIndices.size(); ++trialElemPosition) {
      const int activeTrialElemIndex = trialDeviceElemIndices[trialElemPosition];
      thrust::replace(d_elementPairTrialIndices.begin(), d_elementPairTrialIndices.end(),
          activeTrialElemIndex, trialElemPosition);
    }
  }

  // Allocate device memory for the result
  thrust::device_vector<double> d_result(
      geometryPairCount * testDofCount * trialDofCount);

  // Setup geometry data for selected elements on the device
  m_testGrid->setupElements(testDeviceElemIndices);
  m_trialGrid->setupElements(trialDeviceElemIndices);

  // Measure time of the GPU execution (CUDA event based)
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if (m_cacheElemData == true) {

    thrust::device_vector<double> d_testGeomData;
    thrust::device_vector<double> d_trialGeomData;
    thrust::device_vector<double> d_testNormals;
    thrust::device_vector<double> d_trialNormals;
    thrust::device_vector<double> d_testIntegrationElements;
    thrust::device_vector<double> d_trialIntegrationElements;

    // Calculate global points on the device
    m_testGrid->local2global(
        m_localTestQuadPoints.transpose().eval(), d_testGeomData);
    if (m_testGrid.get() == m_trialGrid.get() &&
        m_localTestQuadPoints == m_localTrialQuadPoints) {
      // TODO: avoid copy
      d_trialGeomData = d_testGeomData;
    } else {
      m_trialGrid->local2global(
          m_localTrialQuadPoints.transpose().eval(), d_trialGeomData);
    }

    m_testGrid->calculateNormalsAndIntegrationElements(
        d_testNormals, d_testIntegrationElements);
    m_trialGrid->calculateNormalsAndIntegrationElements(
        d_trialNormals, d_trialIntegrationElements);

    const unsigned int activeTestElemCount = d_testIntegrationElements.size();
    const unsigned int activeTrialElemCount = d_trialIntegrationElements.size();

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
                     evaluateIntegralFunctor(
                         d_elementPairTestIndices.data(),
                         d_elementPairTrialIndices.data(),
                         testPointCount,
                         trialPointCount,
                         testDofCount,
                         trialDofCount,
                         d_testQuadWeights.data(),
                         d_trialQuadWeights.data(),
                         d_testBasisData.data(),
                         d_trialBasisData.data(),
                         d_testGeomData.data(),
                         d_trialGeomData.data(),
                         activeTestElemCount,
                         activeTrialElemCount,
                         d_testNormals.data(),
                         d_trialNormals.data(),
                         d_testIntegrationElements.data(),
                         d_trialIntegrationElements.data(),
                         d_result.data()
                         ));
  } else {

    thrust::device_vector<double> d_testVtx0x;
    thrust::device_vector<double> d_testVtx0y;
    thrust::device_vector<double> d_testVtx0z;
    thrust::device_vector<double> d_testVtx1x;
    thrust::device_vector<double> d_testVtx1y;
    thrust::device_vector<double> d_testVtx1z;
    thrust::device_vector<double> d_testVtx2x;
    thrust::device_vector<double> d_testVtx2y;
    thrust::device_vector<double> d_testVtx2z;

    thrust::device_vector<double> d_trialVtx0x;
    thrust::device_vector<double> d_trialVtx0y;
    thrust::device_vector<double> d_trialVtx0z;
    thrust::device_vector<double> d_trialVtx1x;
    thrust::device_vector<double> d_trialVtx1y;
    thrust::device_vector<double> d_trialVtx1z;
    thrust::device_vector<double> d_trialVtx2x;
    thrust::device_vector<double> d_trialVtx2y;
    thrust::device_vector<double> d_trialVtx2z;

    m_testGrid->getRawElementData(d_testVtx0x, d_testVtx0y, d_testVtx0z,
                                  d_testVtx1x, d_testVtx1y, d_testVtx1z,
                                  d_testVtx2x, d_testVtx2y, d_testVtx2z);

    m_trialGrid->getRawElementData(d_trialVtx0x, d_trialVtx0y, d_trialVtx0z,
                                   d_trialVtx1x, d_trialVtx1y, d_trialVtx1z,
                                   d_trialVtx2x, d_trialVtx2y, d_trialVtx2z);

    const unsigned int activeTestElemCount = d_testVtx0x.size();
    const unsigned int activeTrialElemCount = d_trialVtx0x.size();

    const unsigned int testPointDim = m_localTestQuadPoints.rows();
    const unsigned int trialPointDim = m_localTrialQuadPoints.rows();

    if (testPointDim != 2 || trialPointDim != 2)
      throw std::runtime_error("CudaIntegrator::integrate(): "
                               "only valid for two-dimensional local points");

    Matrix<double> localTestQuadPointsTransposed = m_localTestQuadPoints.transpose().eval();
    Matrix<double> localTrialQuadPointsTransposed = m_localTrialQuadPoints.transpose().eval();
    thrust::host_vector<double> h_localTestPoints(
        localTestQuadPointsTransposed.data(),
        localTestQuadPointsTransposed.data()+testPointDim*testPointCount);
    thrust::host_vector<double> h_localTrialPoints(
        localTrialQuadPointsTransposed.data(),
        localTrialQuadPointsTransposed.data()+trialPointDim*trialPointCount);

    thrust::host_vector<double> h_testFun0(testPointCount);
    thrust::host_vector<double> h_testFun1(testPointCount);
    thrust::host_vector<double> h_testFun2(testPointCount);
    thrust::host_vector<double> h_trialFun0(trialPointCount);
    thrust::host_vector<double> h_trialFun1(trialPointCount);
    thrust::host_vector<double> h_trialFun2(trialPointCount);

    // Evaluate geometrical shape function values on the host
    thrust::transform(thrust::host,
      thrust::make_zip_iterator(
        thrust::make_tuple(h_localTestPoints.begin(),
                           h_localTestPoints.begin()+testPointCount)),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_localTestPoints.begin()+testPointCount,
                           h_localTestPoints.end())),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_testFun0.begin(),
                           h_testFun1.begin(),
                           h_testFun2.begin())),
      geomShapeFunFunctor());

    thrust::transform(thrust::host,
      thrust::make_zip_iterator(
        thrust::make_tuple(h_localTrialPoints.begin(),
                           h_localTrialPoints.begin()+trialPointCount)),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_localTrialPoints.begin()+trialPointCount,
                           h_localTrialPoints.end())),
      thrust::make_zip_iterator(
        thrust::make_tuple(h_trialFun0.begin(),
                           h_trialFun1.begin(),
                           h_trialFun2.begin())),
      geomShapeFunFunctor());

    // Copy data to device
    thrust::device_vector<double> d_testFun0 = h_testFun0;
    thrust::device_vector<double> d_testFun1 = h_testFun1;
    thrust::device_vector<double> d_testFun2 = h_testFun2;
    thrust::device_vector<double> d_trialFun0 = h_trialFun0;
    thrust::device_vector<double> d_trialFun1 = h_trialFun1;
    thrust::device_vector<double> d_trialFun2 = h_trialFun2;

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
                     evaluateIntegralFunctorNonCached(
                         d_elementPairTestIndices.data(),
                         d_elementPairTrialIndices.data(),
                         testPointCount,
                         trialPointCount,
                         testDofCount,
                         trialDofCount,
                         d_testQuadWeights.data(),
                         d_trialQuadWeights.data(),
                         d_testBasisData.data(),
                         d_trialBasisData.data(),
                         activeTestElemCount,
                         activeTrialElemCount,
                         d_testVtx0x.data(),
                         d_testVtx0y.data(),
                         d_testVtx0z.data(),
                         d_testVtx1x.data(),
                         d_testVtx1y.data(),
                         d_testVtx1z.data(),
                         d_testVtx2x.data(),
                         d_testVtx2y.data(),
                         d_testVtx2z.data(),
                         d_trialVtx0x.data(),
                         d_trialVtx0y.data(),
                         d_trialVtx0z.data(),
                         d_trialVtx1x.data(),
                         d_trialVtx1y.data(),
                         d_trialVtx1z.data(),
                         d_trialVtx2x.data(),
                         d_trialVtx2y.data(),
                         d_trialVtx2z.data(),
                         d_testFun0.data(),
                         d_testFun1.data(),
                         d_testFun2.data(),
                         d_trialFun0.data(),
                         d_trialFun1.data(),
                         d_trialFun2.data(),
                         d_result.data()
                         ));
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeIntegral;
  cudaEventElapsedTime(&elapsedTimeIntegral , start, stop);
  std::cout << "Time for integral evaluation is "
    << elapsedTimeIntegral << " ms" << std::endl;

  cudaEventRecord(start, 0);

  // Copy result back to host memory
  thrust::host_vector<double> h_result = d_result;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeResultCopy;
  cudaEventElapsedTime(&elapsedTimeResultCopy, start, stop);
  std::cout << "Time for result copy is "
    << elapsedTimeResultCopy << " ms" << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Assemble result
  for (int geometryPair = 0; geometryPair < geometryPairCount; ++geometryPair) {
    for (int testDof = 0; testDof < testDofCount; ++testDof) {
      for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
        (*result[geometryPair])(testDof, trialDof) =
            h_result[geometryPair * testDofCount * trialDofCount
                   + testDof * trialDofCount
                   + trialDof];
      }
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for local result assembly = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaIntegrator);

} // namespace Fiber
