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

#include "../common/not_implemented_error.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename ResultType>
CudaIntegrator<BasisFunctionType, ResultType>::CudaIntegrator(
    const Matrix<double> &localTestQuadPoints,
    const Matrix<double> &localTrialQuadPoints,
    const std::vector<double> &testQuadWeights,
    const std::vector<double> &trialQuadWeights,
    shared_ptr<const CudaGrid> testGrid,
    shared_ptr<const CudaGrid> trialGrid)
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
    std::vector<Matrix<ResultType>> &result) const {

  throw Bempp::NotImplementedError(
      "CudaIntegrator::integrate(): not implemented yet");
}

} // namespace Fiber
