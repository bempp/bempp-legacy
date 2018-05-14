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
#include "tbb/enumerable_thread_specific.h"
#include "tbb/atomic.h"

namespace {

  template<typename ResultType, typename CudaResultType>
  struct CudaHostDeviceBuffer {

    CudaHostDeviceBuffer() {

      throw std::invalid_argument(
          "CudaHostDeviceBuffer::CudaHostDeviceBuffer(): "
          "Default constructor called.");
    }

    CudaHostDeviceBuffer(const int deviceId, const size_t maxActiveElemPairCount,
                         const unsigned int trialDofCount, const unsigned int testDofCount) {

      size_t resultArraySize = maxActiveElemPairCount
          * trialDofCount * testDofCount * sizeof(CudaResultType);
      if (Bempp::ScalarTraits<ResultType>::isComplex) resultArraySize *= 2;

      // Allocate pageable host memory and device memory
      cudaSetDevice(deviceId);
      h_resultEven = (CudaResultType*)malloc(resultArraySize);
      h_resultOdd  = (CudaResultType*)malloc(resultArraySize);
      cu_verify( cudaMalloc((void**)&(d_resultEven), resultArraySize) );
      cu_verify( cudaMalloc((void**)&(d_resultOdd ), resultArraySize) );
//      printf("Thread %d allocated memory on device %d.\n", tbb::this_tbb_thread::get_id(), deviceId);
    }

    ~CudaHostDeviceBuffer() {

      // Free pageable host memory and device memory
//      cudaSetDevice(0);
      if (h_resultEven != nullptr) free(h_resultEven);
      if (h_resultOdd  != nullptr) free(h_resultOdd );
      if (d_resultEven != nullptr) cu_verify( cudaFree(d_resultEven) );
      if (d_resultOdd  != nullptr) cu_verify( cudaFree(d_resultOdd ) );
//      printf("Thread %d freed memory on the device.\n", tbb::this_tbb_thread::get_id());
    }

    CudaResultType* h_resultEven;
    CudaResultType* d_resultEven;
    CudaResultType* h_resultOdd;
    CudaResultType* d_resultOdd;
  };
}

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
      const Bempp::Context<BasisFunctionType, ResultType> &context);

  virtual ~CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces();

  // Old implementation with integration on GPU and assembly on CPU
  // Each element pair is treated by one GPU thread.
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
      Matrix<ResultType> &result, const unsigned int device/*,
      tbb::atomic<std::chrono::steady_clock::duration>& allocationTimer,
      tbb::atomic<std::chrono::steady_clock::duration>& integrationTimer,
      tbb::atomic<std::chrono::steady_clock::duration>& assemblyTimer*/);

  // New implementation with both integration and assembly on GPU
  // Each global dof is treated by one GPU thread.
  // - More input data loaded from global device memory
  // - More kernel function evaluations
  // - More input data copied to the device
  // + Avoids computation of unutilized local dofs
  // + Less output data written to global device memory
  // + Avoids slow assembly on the CPU
  // + Less host and device memory allocated
  // + Less output data copied from the device
  virtual void
  evaluateLocalWeakForms(const std::vector<std::vector<Bempp::LocalDof>> & testRawLocalDofs,
                         const std::vector<std::vector<Bempp::LocalDof>> &trialRawLocalDofs,
                         Matrix<ResultType> &result, const unsigned int device);

private:
  /** \cond PRIVATE */
  std::vector<int> m_deviceIds;
  std::vector<shared_ptr<Integrator>> m_cudaIntegrators;
  size_t m_maxActiveElemPairCount;
  const unsigned int m_trialDofCount;
  const unsigned int m_testDofCount;
  std::vector< tbb::enumerable_thread_specific<
      ::CudaHostDeviceBuffer<ResultType, CudaResultType> > > m_buffer;
  /** \endcond */
};

} // namespace Fiber

#include "cuda_default_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

#endif
