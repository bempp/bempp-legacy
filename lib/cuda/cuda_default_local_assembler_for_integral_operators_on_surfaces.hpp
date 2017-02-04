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

#ifndef fiber_cuda_default_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_cuda_default_local_assembler_for_integral_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "cuda_integrator.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
class CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces {

public:

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
  typedef typename ScalarTraits<CudaResultType>::RealType CudaCoordinateType;

  typedef CudaIntegrator<BasisFunctionType, KernelType, ResultType,
      CudaBasisFunctionType, CudaKernelType, CudaResultType> Integrator;

  CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces(
      const Bempp::Space<BasisFunctionType> &testSpace,
      const Bempp::Space<BasisFunctionType> &trialSpace,
      const shared_ptr<const CollectionOfKernels<KernelType>> &kernel,
      const Shapeset<BasisFunctionType> &testShapeset,
      const Shapeset<BasisFunctionType> &trialShapeset,
      const Bempp::CudaOptions &cudaOptions);

  virtual ~CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces();

  virtual void
  evaluateLocalWeakForms(
      const std::vector<int> &testElementIndices,
      const std::vector<int> &trialElementIndices,
      const std::vector<std::vector<LocalDofIndex>> &testLocalDofs,
      const std::vector<std::vector<LocalDofIndex>> &trialLocalDofs,
      const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
      const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
      const std::vector<std::vector<int>> &blockRows,
      const std::vector<std::vector<int>> &blockCols,
      Matrix<ResultType> &result, const unsigned int device = 0);

private:
  /** \cond PRIVATE */
  std::vector<int> m_deviceIds;
  std::vector<shared_ptr<Integrator>> m_cudaIntegrators;
  size_t m_maxActiveElemPairCount;
  const unsigned int m_trialDofCount;
  const unsigned int m_testDofCount;
  /** \endcond */
};

} // namespace Fiber

#include "cuda_default_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

#endif
