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

#ifndef fiber_cuda_integrator_hpp
#define fiber_cuda_integrator_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "../fiber/types.hpp"
#include "../fiber/shared_ptr.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Shapeset;
class CudaGrid;
/** \endcond */

/** \brief Regular integration over pairs of elements on the device. */
template <typename BasisFunctionType, typename ResultType>
class CudaIntegrator {
public:

  /** \brief Constructor */
  CudaIntegrator(
      const Matrix<double> &localTestQuadPoints,
      const Matrix<double> &localTrialQuadPoints,
      const std::vector<double> &testQuadWeights,
      const std::vector<double> &trialQuadWeights,
      shared_ptr<const CudaGrid> testGrid,
      shared_ptr<const CudaGrid> trialGrid);

  /** \brief Destructor. */
  virtual ~CudaIntegrator();

  void integrate(
      const std::vector<int> &elementPairTestIndices,
      const std::vector<int> &elementPairTrialIndices,
      const Shapeset<BasisFunctionType> &testShapeset,
      const Shapeset<BasisFunctionType> &trialShapeset,
      std::vector<Matrix<ResultType>> &result) const;

private:
  /** \cond PRIVATE */

  Matrix<double> m_localTestQuadPoints;
  Matrix<double> m_localTrialQuadPoints;
  std::vector<double> m_testQuadWeights;
  std::vector<double> m_trialQuadWeights;

  shared_ptr<const CudaGrid> m_testGrid;
  shared_ptr<const CudaGrid> m_trialGrid;

  /** \endcond */
};

} // namespace Fiber

#endif
