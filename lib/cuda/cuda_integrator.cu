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

struct BasisFunData {
  unsigned int dofCount;
  thrust::device_ptr<double> basisData;
};

struct QuadData {
  unsigned int pointCount;
  thrust::device_ptr<double> weights;
};

struct GeomShapeFunData {
  thrust::device_ptr<const double> fun0;
  thrust::device_ptr<const double> fun1;
  thrust::device_ptr<const double> fun2;
};

struct RawGeometryData {
  unsigned int elemCount, vtxCount;
  thrust::device_ptr<const double> vertices;
  thrust::device_ptr<const int> elementCorners;
};

struct ElemData {
  unsigned int activeElemCount;
  thrust::device_ptr<double> geomData;
  thrust::device_ptr<const double> normals;
  thrust::device_ptr<const double> integrationElements;
};

//template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct evaluateIntegralFunctorCached {

  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  QuadData testQuadData, trialQuadData;
  unsigned int testDofCount, trialDofCount;
  thrust::device_ptr<double> testBasisData, trialBasisData;
  ElemData testElemData, trialElemData;
  thrust::device_ptr<double> result;

  evaluateIntegralFunctorCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData _testQuadData, const QuadData _trialQuadData,
      const unsigned int _testDofCount, const unsigned int _trialDofCount,
      thrust::device_ptr<double> _testBasisData,
      thrust::device_ptr<double> _trialBasisData,
      const ElemData _testElemData, const ElemData _trialElemData,
      thrust::device_ptr<double> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testQuadData(_testQuadData), trialQuadData(_trialQuadData),
        testDofCount(_testDofCount), trialDofCount(_trialDofCount),
        testBasisData(_testBasisData), trialBasisData(_trialBasisData),
        testElemData(_testElemData), trialElemData(_trialElemData),
        result (_result) {

    assert(testQuadData.pointCount <= 10);
    assert(trialQuadData.pointCount <= 10);
    assert(testDofCount * trialDofCount <= 36);
  }

  __host__ __device__
  void operator() (const unsigned int i) {

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    clock_t start_time = clock();

    const double testIntegrationElement =
        testElemData.integrationElements[testElemPosition];
    const double trialIntegrationElement =
        trialElemData.integrationElements[trialElemPosition];

    const unsigned int activeTestElemCount = testElemData.activeElemCount;
    const unsigned int activeTrialElemCount = trialElemData.activeElemCount;

    const double xTestElemNormal = testElemData.normals[testElemPosition];
    const double yTestElemNormal = testElemData.normals[testElemPosition+activeTestElemCount];
    const double zTestElemNormal = testElemData.normals[testElemPosition+2*activeTestElemCount];

    const double xTrialElemNormal = trialElemData.normals[trialElemPosition];
    const double yTrialElemNormal = trialElemData.normals[trialElemPosition+activeTrialElemCount];
    const double zTrialElemNormal = trialElemData.normals[trialElemPosition+2*activeTrialElemCount];

    clock_t stop_time = clock();

    if (i == 0)
      printf("Gather normals and integration elements: %d clock cycles\n",
          stop_time-start_time);

    start_time = clock();

    double xTestGeomData[10], yTestGeomData[10], zTestGeomData[10];
    double xTrialGeomData[10], yTrialGeomData[10], zTrialGeomData[10];
    double testElemQuadWeights[10], trialElemQuadWeights[10];

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
            const double xDist = xTestGeomData[testPoint] - xTrialGeomData[trialPoint];
            const double yDist = yTestGeomData[testPoint] - yTrialGeomData[trialPoint];
            const double zDist = zTestGeomData[testPoint] - zTrialGeomData[trialPoint];
            const double distance = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
            const double kernelValue = 1.0 / (4.0 * M_PI * distance);
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
  }
};

//template <typename BasisFunctionType, typename KernelType, typename ResultType>
struct evaluateIntegralFunctorNonCached {

  thrust::device_ptr<int> elementPairTestIndexPositions;
  thrust::device_ptr<int> elementPairTrialIndexPositions;
  QuadData testQuadData, trialQuadData;
  unsigned int testDofCount, trialDofCount;
  thrust::device_ptr<double> testBasisData, trialBasisData;
  RawGeometryData testRawGeometryData, trialRawGeometryData;
  GeomShapeFunData testGeomShapeFunData, trialGeomShapeFunData;
  thrust::device_ptr<double> result;

  evaluateIntegralFunctorNonCached(
      thrust::device_ptr<int> _elementPairTestIndexPositions,
      thrust::device_ptr<int> _elementPairTrialIndexPositions,
      const QuadData _testQuadData, const QuadData _trialQuadData,
      const unsigned int _testDofCount, const unsigned int _trialDofCount,
      thrust::device_ptr<double> _testBasisData,
      thrust::device_ptr<double> _trialBasisData,
      const RawGeometryData _testRawGeometryData,
      const RawGeometryData _trialRawGeometryData,
      const GeomShapeFunData _testGeomShapeFunData,
      const GeomShapeFunData _trialGeomShapeFunData,
      thrust::device_ptr<double> _result)
      : elementPairTestIndexPositions(_elementPairTestIndexPositions),
        elementPairTrialIndexPositions(_elementPairTrialIndexPositions),
        testQuadData(_testQuadData), trialQuadData(_trialQuadData),
        testDofCount(_testDofCount), trialDofCount(_trialDofCount),
        testBasisData(_testBasisData), trialBasisData(_trialBasisData),
        testRawGeometryData(_testRawGeometryData),
        trialRawGeometryData(_trialRawGeometryData),
        testGeomShapeFunData(_testGeomShapeFunData),
        trialGeomShapeFunData(_trialGeomShapeFunData),
        result (_result) {

    assert(testQuadData.pointCount <= 10);
    assert(trialQuadData.pointCount <= 10);
    assert(testDofCount * trialDofCount <= 36);
  }

  __host__ __device__
  void operator() (const unsigned int i) {

    const int coordCount = 3;

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    clock_t start_time = clock();

    double testElemVtx0[coordCount];
    double testElemVtx1[coordCount];
    double testElemVtx2[coordCount];

    double trialElemVtx0[coordCount];
    double trialElemVtx1[coordCount];
    double trialElemVtx2[coordCount];

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

    double xTestElemNormal =
        (testElemVtx1[1] - testElemVtx0[1]) * (testElemVtx2[2] - testElemVtx0[2])
      - (testElemVtx1[2] - testElemVtx0[2]) * (testElemVtx2[1] - testElemVtx0[1]);
    double yTestElemNormal =
        (testElemVtx1[2] - testElemVtx0[2]) * (testElemVtx2[0] - testElemVtx0[0])
      - (testElemVtx1[0] - testElemVtx0[0]) * (testElemVtx2[2] - testElemVtx0[2]);
    double zTestElemNormal =
        (testElemVtx1[0] - testElemVtx0[0]) * (testElemVtx2[1] - testElemVtx0[1])
      - (testElemVtx1[1] - testElemVtx0[1]) * (testElemVtx2[0] - testElemVtx0[0]);

    double xTrialElemNormal =
        (trialElemVtx1[1] - trialElemVtx0[1]) * (trialElemVtx2[2] - trialElemVtx0[2])
      - (trialElemVtx1[2] - trialElemVtx0[2]) * (trialElemVtx2[1] - trialElemVtx0[1]);
    double yTrialElemNormal =
        (trialElemVtx1[2] - trialElemVtx0[2]) * (trialElemVtx2[0] - trialElemVtx0[0])
      - (trialElemVtx1[0] - trialElemVtx0[0]) * (trialElemVtx2[2] - trialElemVtx0[2]);
    double zTrialElemNormal =
        (trialElemVtx1[0] - trialElemVtx0[0]) * (trialElemVtx2[1] - trialElemVtx0[1])
      - (trialElemVtx1[1] - trialElemVtx0[1]) * (trialElemVtx2[0] - trialElemVtx0[0]);

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

    double testGeomData[10 * coordCount], trialGeomData[10 * coordCount];
    double testElemQuadWeights[10], trialElemQuadWeights[10];

    const unsigned int testPointCount = testQuadData.pointCount;
    const unsigned int trialPointCount = trialQuadData.pointCount;

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const double ptFun0 = testGeomShapeFunData.fun0[testPoint];
      const double ptFun1 = testGeomShapeFunData.fun1[testPoint];
      const double ptFun2 = testGeomShapeFunData.fun2[testPoint];
      for (int i = 0; i < coordCount; ++i) {
        testGeomData[coordCount * testPoint + i] =
            ptFun0 * testElemVtx0[i]
          + ptFun1 * testElemVtx1[i]
          + ptFun2 * testElemVtx2[i];
      }
      testElemQuadWeights[testPoint] = testQuadData.weights[testPoint];
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const double ptFun0 = trialGeomShapeFunData.fun0[trialPoint];
      const double ptFun1 = trialGeomShapeFunData.fun1[trialPoint];
      const double ptFun2 = trialGeomShapeFunData.fun2[trialPoint];
      for (int i = 0; i < coordCount; ++i) {
        trialGeomData[coordCount * trialPoint + i] =
            ptFun0 * trialElemVtx0[i]
          + ptFun1 * trialElemVtx1[i]
          + ptFun2 * trialElemVtx2[i];
      }
      trialElemQuadWeights[trialPoint] = trialQuadData.weights[trialPoint];
    }

    stop_time = clock();

    if (i == 0)
      printf("Calculate global points: %d clock cycles\n", stop_time-start_time);

    start_time = clock();

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
            const double xDist = testGeomData[coordCount * testPoint]
                               - trialGeomData[coordCount * trialPoint];
            const double yDist = testGeomData[coordCount * testPoint + 1]
                               - trialGeomData[coordCount * trialPoint + 1];
            const double zDist = testGeomData[coordCount * testPoint + 2]
                               - trialGeomData[coordCount * trialPoint + 2];
            const double distance = std::sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
            const double kernelValue = 1.0 / (4.0 * M_PI * distance);
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
    const Shapeset<BasisFunctionType> &testShapeset,
    const Shapeset<BasisFunctionType> &trialShapeset,
    std::vector<Matrix<ResultType>*> &result) {

  std::cout << "Hello, this is CudaIntegrator::integrate()!" << std::endl;

  const unsigned int testPointCount = m_localTestQuadPoints.cols();
  const unsigned int trialPointCount = m_localTrialQuadPoints.cols();

  const unsigned int testPointDim = m_localTestQuadPoints.rows();
  const unsigned int trialPointDim = m_localTrialQuadPoints.rows();

  if (testPointDim != 2 || trialPointDim != 2)
    throw std::invalid_argument("CudaIntegrator::integrate(): "
                                "only valid for two-dimensional local points");

  const unsigned int geometryPairCount = elementPairTestIndices.size();

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
  QuadData testQuadData, trialQuadData;
  thrust::device_vector<double> d_testQuadWeights(m_testQuadWeights);
  thrust::device_vector<double> d_trialQuadWeights;
  testQuadData.weights = d_testQuadWeights.data();
  if (m_testQuadWeights == m_trialQuadWeights) {
    trialQuadData.weights = testQuadData.weights;
  } else {
    d_trialQuadWeights.resize(m_trialQuadWeights.size());
    d_trialQuadWeights.assign(
        m_trialQuadWeights.begin(), m_trialQuadWeights.end());
    trialQuadData.weights = d_trialQuadWeights.data();
  }
  testQuadData.pointCount = testPointCount;
  trialQuadData.pointCount = trialPointCount;

  // Copy element pair indices to device memory
  thrust::device_vector<int> d_elementPairTestIndices(elementPairTestIndices);
  thrust::device_vector<int> d_elementPairTrialIndices(elementPairTrialIndices);

  // Allocate device memory for the result
  thrust::device_vector<double> d_result(
      geometryPairCount * testDofCount * trialDofCount);

  // Measure time of the GPU execution (CUDA event based)
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  if (m_cacheElemData == true) {

    // Setup geometry data for selected elements on the device
    m_testGrid->setupGeometry();
    m_trialGrid->setupGeometry();

    ElemData testElemData, trialElemData;

    thrust::device_vector<double> d_testGeomData, d_trialGeomData;
    thrust::device_vector<double> d_testNormals, d_trialNormals;
    thrust::device_vector<double> d_testIntegrationElements, d_trialIntegrationElements;

    // Calculate global points on the device
    m_testGrid->local2global(m_localTestQuadPoints, d_testGeomData);
    testElemData.geomData = d_testGeomData.data();
    if (m_testGrid.get() == m_trialGrid.get() &&
        m_localTestQuadPoints == m_localTrialQuadPoints) {
      trialElemData.geomData = testElemData.geomData;
    } else {
      m_trialGrid->local2global(m_localTrialQuadPoints, d_trialGeomData);
      trialElemData.geomData = d_trialGeomData.data();
    }

    m_testGrid->calculateNormalsAndIntegrationElements(
        d_testNormals, d_testIntegrationElements);
    testElemData.normals = d_testNormals.data();
    testElemData.integrationElements = d_testIntegrationElements.data();
    testElemData.activeElemCount = d_testIntegrationElements.size();
    if (m_testGrid.get() == m_trialGrid.get()) {
      trialElemData.normals = testElemData.normals;
      trialElemData.integrationElements = testElemData.integrationElements;
      trialElemData.activeElemCount = testElemData.activeElemCount;
    } else {
      m_trialGrid->calculateNormalsAndIntegrationElements(
          d_trialNormals, d_trialIntegrationElements);
      trialElemData.normals = d_trialNormals.data();
      trialElemData.integrationElements = d_trialIntegrationElements.data();
      trialElemData.activeElemCount = d_trialIntegrationElements.size();
    }

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
                     evaluateIntegralFunctorCached(
                         d_elementPairTestIndices.data(),
                         d_elementPairTrialIndices.data(),
                         testQuadData, trialQuadData,
                         testDofCount, trialDofCount,
                         d_testBasisData.data(), d_trialBasisData.data(),
                         testElemData, trialElemData,
                         d_result.data()
                         ));
  } else {

    // Get raw geometry data
    RawGeometryData testRawGeometryData, trialRawGeometryData;
    m_testGrid->getRawGeometryData(
        testRawGeometryData.vtxCount, testRawGeometryData.elemCount,
        testRawGeometryData.vertices, testRawGeometryData.elementCorners);
    m_trialGrid->getRawGeometryData(
        trialRawGeometryData.vtxCount, trialRawGeometryData.elemCount,
        trialRawGeometryData.vertices, trialRawGeometryData.elementCorners);

    thrust::host_vector<double> h_testFun0(testPointCount);
    thrust::host_vector<double> h_testFun1(testPointCount);
    thrust::host_vector<double> h_testFun2(testPointCount);
    thrust::host_vector<double> h_trialFun0(trialPointCount);
    thrust::host_vector<double> h_trialFun1(trialPointCount);
    thrust::host_vector<double> h_trialFun2(trialPointCount);

    // Evaluate geometrical shape function values
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const double r = m_localTestQuadPoints(0,testPoint);
      const double s = m_localTestQuadPoints(1,testPoint);
      h_testFun0[testPoint] = 1.0 - r - s;
      h_testFun1[testPoint] = r;
      h_testFun2[testPoint] = s;
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const double r = m_localTrialQuadPoints(0,trialPoint);
      const double s = m_localTrialQuadPoints(1,trialPoint);
      h_trialFun0[trialPoint] = 1.0 - r - s;
      h_trialFun1[trialPoint] = r;
      h_trialFun2[trialPoint] = s;
    }

    // Copy data to device
    thrust::device_vector<double> d_testFun0 = h_testFun0;
    thrust::device_vector<double> d_testFun1 = h_testFun1;
    thrust::device_vector<double> d_testFun2 = h_testFun2;
    thrust::device_vector<double> d_trialFun0 = h_trialFun0;
    thrust::device_vector<double> d_trialFun1 = h_trialFun1;
    thrust::device_vector<double> d_trialFun2 = h_trialFun2;

    // Collect geometrical shape function data
    GeomShapeFunData testGeomShapeFunData, trialGeomShapeFunData;
    testGeomShapeFunData.fun0 = d_testFun0.data();
    testGeomShapeFunData.fun1 = d_testFun1.data();
    testGeomShapeFunData.fun2 = d_testFun2.data();
    trialGeomShapeFunData.fun0 = d_trialFun0.data();
    trialGeomShapeFunData.fun1 = d_trialFun1.data();
    trialGeomShapeFunData.fun2 = d_trialFun2.data();

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
                     evaluateIntegralFunctorNonCached(
                         d_elementPairTestIndices.data(),
                         d_elementPairTrialIndices.data(),
                         testQuadData, trialQuadData,
                         testDofCount, trialDofCount,
                         d_testBasisData.data(), d_trialBasisData.data(),
                         testRawGeometryData, trialRawGeometryData,
                         testGeomShapeFunData, trialGeomShapeFunData,
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

//  std::cout << "h_result = " << std::endl;
//  for (int i = 0; i < h_result.size(); ++i) {
//    std::cout << h_result[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

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
