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
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>

namespace Fiber {

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
        result (_result) { }

  __host__ __device__
  void operator() (const unsigned int i) {

    const int testElemPosition = elementPairTestIndexPositions[i];
    const int trialElemPosition = elementPairTrialIndexPositions[i];

    const double testIntegrationElement = testIntegrationElements[testElemPosition];
    const double trialIntegrationElement = trialIntegrationElements[trialElemPosition];

    const double xTestElemNormal = testNormals[testElemPosition];
    const double yTestElemNormal = testNormals[testElemPosition+activeTestElemCount];
    const double zTestElemNormal = testNormals[testElemPosition+2*activeTestElemCount];

    const double xTrialElemNormal = trialNormals[trialElemPosition];
    const double yTrialElemNormal = trialNormals[trialElemPosition+activeTrialElemCount];
    const double zTrialElemNormal = trialNormals[trialElemPosition+2*activeTrialElemCount];

    double* xTestGeomData = new double[testPointCount];
    double* yTestGeomData = new double[testPointCount];
    double* zTestGeomData = new double[testPointCount];

    double* xTrialGeomData = new double[trialPointCount];
    double* yTrialGeomData = new double[trialPointCount];
    double* zTrialGeomData = new double[trialPointCount];

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      xTestGeomData[testPoint] =
          testGeomData[testPoint * activeTestElemCount + testElemPosition];
      yTestGeomData[testPoint] =
          testGeomData[testPoint * activeTestElemCount + testElemPosition
                       + testPointCount * activeTestElemCount];
      zTestGeomData[testPoint] =
          testGeomData[testPoint * activeTestElemCount + testElemPosition
                       + 2 * testPointCount * activeTestElemCount];
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
    }

    double* localResult = new double[testDofCount * trialDofCount];
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        double sum = 0.;
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
          const double trialWeight =
              trialIntegrationElement * trialQuadWeights[trialPoint];
          double partialSum = 0.;
          for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
            const double testWeight =
                testIntegrationElement * testQuadWeights[testPoint];
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

    // Copy local result to global device memory
    const unsigned int offset = i * testDofCount * trialDofCount;
    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
        result[offset + testDof * trialDofCount + trialDof] =
            localResult[testDof * trialDofCount + trialDof];
      }
    }

    // Free local arrays
    delete[] xTestGeomData;
    delete[] yTestGeomData;
    delete[] zTestGeomData;

    delete[] xTrialGeomData;
    delete[] yTrialGeomData;
    delete[] zTrialGeomData;
  }
};

template <typename BasisFunctionType, typename ResultType>
CudaIntegrator<BasisFunctionType, ResultType>::CudaIntegrator(
    const Matrix<double> &localTestQuadPoints,
    const Matrix<double> &localTrialQuadPoints,
    const std::vector<double> &testQuadWeights,
    const std::vector<double> &trialQuadWeights,
    shared_ptr<const Bempp::CudaGrid> testGrid,
    shared_ptr<const Bempp::CudaGrid> trialGrid)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testQuadWeights(testQuadWeights), m_trialQuadWeights(trialQuadWeights),
      m_testGrid(testGrid), m_trialGrid(trialGrid) {

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
    std::vector<Matrix<ResultType>*> &result) const {

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

  // Calculate global points on the device
  thrust::device_vector<double> d_testGeomData;
  thrust::device_vector<double> d_trialGeomData;
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

  unsigned int activeTestElemCount;
  unsigned int activeTrialElemCount;
  thrust::device_vector<int> d_activeTestElemIndices;
  thrust::device_vector<int> d_activeTrialElemIndices;
  thrust::device_vector<double> d_testNormals;
  thrust::device_vector<double> d_trialNormals;
  thrust::device_vector<double> d_testIntegrationElements;
  thrust::device_vector<double> d_trialIntegrationElements;
  m_testGrid->getElementData(activeTestElemCount, d_activeTestElemIndices,
      d_testNormals, d_testIntegrationElements);
  m_trialGrid->getElementData(activeTrialElemCount, d_activeTrialElemIndices,
      d_trialNormals, d_trialIntegrationElements);

  // Copy numerical quadrature weights to device memory
  thrust::device_vector<double> d_testQuadWeights(m_testQuadWeights);
  thrust::device_vector<double> d_trialQuadWeights(m_trialQuadWeights);

  // Copy element pair indices to device memory
  thrust::device_vector<int> d_elementPairTestIndices(elementPairTestIndices);
  thrust::device_vector<int> d_elementPairTrialIndices(elementPairTrialIndices);

  // TODO Replace element pair indices by their positions in activeElemIndices
  // vectors (not necessary in case of ALL_ELEMS active)
  for (int testElemPosition = 0; testElemPosition < activeTestElemCount; ++testElemPosition) {
    const int activeTestElemIndex = d_activeTestElemIndices[testElemPosition];
    thrust::replace(d_elementPairTestIndices.begin(), d_elementPairTestIndices.end(),
        activeTestElemIndex, testElemPosition);
  }
  for (int trialElemPosition = 0; trialElemPosition < activeTrialElemCount; ++trialElemPosition) {
    const int activeTrialElemIndex = d_activeTrialElemIndices[trialElemPosition];
    thrust::replace(d_elementPairTrialIndices.begin(), d_elementPairTrialIndices.end(),
        activeTrialElemIndex, trialElemPosition);
  }

  // Allocate device memory for the result
  thrust::device_vector<double> d_result(
      geometryPairCount * testDofCount * trialDofCount);

  // Measure time of the GPU execution (CUDA event based)
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

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

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeIntegral;
  cudaEventElapsedTime(&elapsedTimeIntegral , start, stop);
  std::cout << "Time for integral evaluation is "
    << elapsedTimeIntegral << " ms" << std::endl;

  // Copy result back to host memory
  thrust::host_vector<double> h_result = d_result;

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
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaIntegrator);

} // namespace Fiber
