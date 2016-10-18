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

#include <thrust/device_vector.h>

namespace Fiber {

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

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  size_t testBasisDeps = 0, trialBasisDeps = 0;
//  testShapeset.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
//                        testBasisData);
//  trialShapeset.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
//                         trialBasisData);

  thrust::device_vector<double> testGeomData;
  thrust::device_vector<double> trialGeomData;
  m_testGrid->local2global(
      m_localTestQuadPoints.transpose().eval(), testGeomData);
  if (m_testGrid.get() == m_trialGrid.get() &&
      m_localTestQuadPoints == m_localTrialQuadPoints) {
    trialGeomData.data() = testGeomData.data();
  } else {
    m_trialGrid->local2global(
        m_localTrialQuadPoints.transpose().eval(), trialGeomData);
  }

  throw Bempp::NotImplementedError(
      "CudaIntegrator::integrate(): not completed yet");
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaIntegrator);

} // namespace Fiber
