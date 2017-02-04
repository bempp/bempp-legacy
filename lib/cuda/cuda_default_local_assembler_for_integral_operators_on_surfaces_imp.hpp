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

#include "cuda_default_local_assembler_for_integral_operators_on_surfaces.hpp"

#include "../fiber/shapeset.hpp"
#include "../fiber/numerical_quadrature.hpp"

#include "../grid/grid.hpp"

#include "../space/space.hpp"

#include <tbb/task_group.h>

#include <chrono>

namespace Fiber {

// Helper functions and classes
namespace {

template <typename CudaResultType, typename ResultType>
struct AssemblyHelper {
  static ResultType value(
      const CudaResultType *h_result,
      const size_t index, const size_t offset) {
    throw std::invalid_argument(
        "CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces::evaluateLocalWeakForms(): "
        "invalid ResultType");
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, float> {
  static float value(
      const CudaResultType *h_result,
      const size_t index, const size_t offset) {
    return static_cast<float>(h_result[index]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, double> {
  static double value(
      const CudaResultType *h_result,
      const size_t index, const size_t offset) {
    return static_cast<double>(h_result[index]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, std::complex<float>> {
  static std::complex<float> value(
      const CudaResultType *h_result,
      const size_t index, const size_t offset) {
    return std::complex<float>(
        h_result[index], h_result[index + offset]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, std::complex<double>> {
  static std::complex<double> value(
      const CudaResultType *h_result,
      const size_t index, const size_t offset) {
    return std::complex<double>(
        h_result[index], h_result[index + offset]);
  }
};

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void IntegrationTask(tbb::task_group &taskGroupDevice,
    const size_t chunk, const bool isLastChunk,
    const size_t elemPairCount, const size_t maxActiveElemPairCount,
    const thrust::device_vector<int> &d_testElementIndices,
    const thrust::device_vector<int> &d_trialElementIndices,
    shared_ptr<Fiber::CudaIntegrator<BasisFunctionType, KernelType, ResultType,
        CudaBasisFunctionType, CudaKernelType, CudaResultType>> &cudaIntegrator,
    CudaResultType *h_result, CudaResultType *d_result, size_t resultArraySize,
    float &integrationTimer) {

  taskGroupDevice.run_and_wait([
       chunk, isLastChunk, elemPairCount, maxActiveElemPairCount,
       &d_testElementIndices, &d_trialElementIndices,
       &cudaIntegrator, h_result, d_result, resultArraySize, &integrationTimer] {

    const size_t elemPairIndexBegin = chunk * maxActiveElemPairCount;

    size_t elemPairIndexEnd;
    if (isLastChunk)
      elemPairIndexEnd = elemPairCount;
    else
      elemPairIndexEnd = (chunk + 1) * maxActiveElemPairCount;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Evaluate regular integrals over selected element pairs
    cudaIntegrator->integrate(d_testElementIndices, d_trialElementIndices,
        elemPairIndexBegin, elemPairIndexEnd, d_result);

    // Copy results to host
    cu_verify( cudaMemcpy(h_result, d_result, resultArraySize,
        cudaMemcpyDeviceToHost) );

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (chunk != 0 && !isLastChunk)
      integrationTimer += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    if (false)
      std::cout << "Time for CudaIntegrator::integrate() = "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                << " ms" << std::endl;
  });
}

template <typename BasisFunctionType, typename ResultType, typename CudaResultType>
void AssemblyTask(tbb::task_group &taskGroupDevice,
    const size_t chunk, const bool isLastChunk, const size_t chunkElemPairCount,
    const size_t maxActiveElemPairCount,
    const unsigned int testDofCount, const unsigned int trialDofCount,
    const std::vector<int> &testElementIndices,
    const std::vector<int> &trialElementIndices,
    const std::vector<std::vector<LocalDofIndex>> &testLocalDofs,
    const std::vector<std::vector<LocalDofIndex>> &trialLocalDofs,
    const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
    const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
    const std::vector<std::vector<int>> &blockRows,
    const std::vector<std::vector<int>> &blockCols,
    CudaResultType *h_result, Matrix<ResultType> &result,
    Matrix<tbb::spin_mutex> &mutex, float &assemblyTimer) {

  typedef tbb::spin_mutex MutexType;

  taskGroupDevice.run([
      chunk, isLastChunk, chunkElemPairCount, maxActiveElemPairCount,
      testDofCount, trialDofCount, &testElementIndices, &trialElementIndices,
      &testLocalDofs, &trialLocalDofs, &testLocalDofWeights, &trialLocalDofWeights,
      &blockRows, &blockCols, h_result, &result, &mutex, &assemblyTimer] {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Parallel global assembly
    tbb::parallel_for(size_t(0), size_t(chunkElemPairCount), [
        chunk, isLastChunk, chunkElemPairCount, maxActiveElemPairCount,
        testDofCount, trialDofCount, &testElementIndices, &trialElementIndices,
        &testLocalDofs, &trialLocalDofs, &testLocalDofWeights, &trialLocalDofWeights,
        &blockRows, &blockCols, h_result, &result, &mutex](size_t chunkElemPair) {

//    for (size_t chunkElemPair = 0; chunkElemPair < chunkElemPairCount; ++chunkElemPair) {

      const size_t offsetResultImag =
          trialDofCount * testDofCount * chunkElemPairCount;

      const size_t offset = chunk * maxActiveElemPairCount;

      const size_t testIndexCount = testElementIndices.size();

      const int trialIndex = (offset + chunkElemPair) / testIndexCount;
      const int testIndex = (offset + chunkElemPair) % testIndexCount;

      // Add the integrals to appropriate entries in the result array
      for (int trialDof = 0; trialDof < trialLocalDofs[trialIndex].size(); ++trialDof) {
        for (int testDof = 0; testDof < testLocalDofs[testIndex].size(); ++testDof) {
          const size_t index = trialLocalDofs[trialIndex][trialDof] * testDofCount * chunkElemPairCount
                               + testLocalDofs[testIndex][testDof] * chunkElemPairCount
                               + chunkElemPair;
          MutexType::scoped_lock lock(mutex(blockRows[testIndex][testDof],
                                            blockCols[trialIndex][trialDof]));
          result(blockRows[testIndex][testDof],
                 blockCols[trialIndex][trialDof]) +=
              conjugate(testLocalDofWeights[testIndex][testDof]) *
              trialLocalDofWeights[trialIndex][trialDof] *
              AssemblyHelper<CudaResultType, ResultType>::
              value(h_result, index, offsetResultImag);
        }
      }
//    }
    });
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (chunk != 0 && !isLastChunk)
      assemblyTimer += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    if (false)
      std::cout << "Time for global result assembly = "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                << " ms" << std::endl;
  });
}

} // namespace

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces(
        const Bempp::Space<BasisFunctionType> &testSpace,
        const Bempp::Space<BasisFunctionType> &trialSpace,
        const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
        const Fiber::Shapeset<BasisFunctionType> &testShapeset,
        const Fiber::Shapeset<BasisFunctionType> &trialShapeset,
        const Bempp::CudaOptions &cudaOptions)
        : m_testDofCount(testShapeset.size()),
          m_trialDofCount(trialShapeset.size()) {

  std::cout << "NOTE: Shapesets have to be identical within one space" << std::endl;
  std::cout << "trialDofCount = " << m_trialDofCount << ", " << std::flush;
  std::cout << "testDofCount = " << m_testDofCount << std::endl;

  // Get numerical quadrature points and weights
  const int trialQuadOrder = cudaOptions.quadOrder();
  const int testQuadOrder = cudaOptions.quadOrder();
  Matrix<CoordinateType> localTrialQuadPoints, localTestQuadPoints;
  std::vector<CoordinateType> trialQuadWeights, testQuadWeights;
  fillSingleQuadraturePointsAndWeights(3, trialQuadOrder,
      localTrialQuadPoints, trialQuadWeights);
  fillSingleQuadraturePointsAndWeights(3, testQuadOrder,
      localTestQuadPoints, testQuadWeights);
  //  std::cout << "trialQuadOrder = " << trialQuadOrder << ", " << std::flush;
  //  std::cout << "testQuadOrder = " << testQuadOrder << std::endl;

  m_deviceIds = cudaOptions.devices();
  const unsigned int deviceCount = m_deviceIds.size();

  m_cudaIntegrators.resize(deviceCount);

  // Get maximum number of element pairs which can be treated on the device
  // simultaneously
  m_maxActiveElemPairCount = cudaOptions.chunkElemPairCount();

  // Let the chunk size be a multiple of the warp size
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  const int warpSize = prop.warpSize;
  m_maxActiveElemPairCount =
      std::max((m_maxActiveElemPairCount / warpSize) * warpSize, size_t(warpSize));
  std::cout << "maxActiveElemPairCount = " << m_maxActiveElemPairCount << std::endl;

  // Loop over all devices
  for (int device = 0; device < deviceCount; ++device) {

    const int deviceId = m_deviceIds[device];
    cu_verify( cudaSetDevice(deviceId) );

    // Push raw grid data to the device
    shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> trialGrid =
        trialSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceId);
    shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> testGrid =
        testSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceId);

    // Create CUDA integrator
    m_cudaIntegrators[device] = boost::make_shared<Integrator>(
        localTestQuadPoints, localTrialQuadPoints,
        testQuadWeights, trialQuadWeights,
        testShapeset, trialShapeset, testGrid, trialGrid,
        m_maxActiveElemPairCount, kernel, deviceId, cudaOptions);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    ~CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces() {

}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    evaluateLocalWeakForms(
        const std::vector<int> &testElementIndices,
        const std::vector<int> &trialElementIndices,
        const std::vector<std::vector<LocalDofIndex>> &testLocalDofs,
        const std::vector<std::vector<LocalDofIndex>> &trialLocalDofs,
        const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
        const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
        const std::vector<std::vector<int>> &blockRows,
        const std::vector<std::vector<int>> &blockCols,
        Matrix<ResultType> &result, const unsigned int device) {

  const int deviceId = m_deviceIds[device];
  cu_verify( cudaSetDevice(deviceId) );

  // Get the number of element pairs to be treated on the specified device
  const size_t elemPairCount = trialElementIndices.size() * testElementIndices.size();
//  std::cout << "elemPairCount = " << elemPairCount << std::endl;

  // Get number of element pair chunks
  const unsigned int chunkCount = static_cast<unsigned int>(
      (elemPairCount-1) / m_maxActiveElemPairCount + 1);

  CudaResultType *h_resultEven, *h_resultOdd, *d_resultEven, *d_resultOdd;
  size_t resultArraySize;
  if (chunkCount == 1) {

    resultArraySize = elemPairCount
        * m_trialDofCount * m_testDofCount * sizeof(CudaResultType);
    if (Bempp::ScalarTraits<ResultType>::isComplex) resultArraySize *= 2;

    // Allocate pageable host memory and memory on the device
    h_resultEven = (CudaResultType*)malloc(resultArraySize);            // Host
    cu_verify( cudaMalloc((void**)&(d_resultEven), resultArraySize) );  // Device

  } else {

    resultArraySize = m_maxActiveElemPairCount
        * m_trialDofCount * m_testDofCount * sizeof(CudaResultType);
    if (Bempp::ScalarTraits<ResultType>::isComplex) resultArraySize *= 2;

    // Allocate pageable host memory and memory on the device
    h_resultEven = (CudaResultType*)malloc(resultArraySize);            // Host
    h_resultOdd = (CudaResultType*)malloc(resultArraySize);
    cu_verify( cudaMalloc((void**)&(d_resultEven), resultArraySize) );  // Device
    cu_verify( cudaMalloc((void**)&(d_resultOdd), resultArraySize) );
  }

  // Copy element indices to device memory
  thrust::device_vector<int> d_testElementIndices, d_trialElementIndices;
  m_cudaIntegrators[device]->pushElementIndicesToDevice(
      testElementIndices, trialElementIndices,
      d_testElementIndices, d_trialElementIndices);

  // Create a mutex matrix for parallel assembly
  typedef tbb::spin_mutex MutexType;
  Matrix<MutexType> mutex(result.rows(), result.cols());

  tbb::task_group taskGroupDevice;
  float integrationTimer = 0., assemblyTimer = 0.;

  // Loop over chunks of element pairs
  for (size_t chunk = 0; chunk < chunkCount; ++chunk) {

//    std::cout << "chunk = " << chunk << std::endl;

    const bool isLastChunk = (chunk == chunkCount - 1);

    CudaResultType *h_result, *d_result;
    if (chunk % 2 == 0) {
      h_result = h_resultEven;
      d_result = d_resultEven;
    } else {
      h_result = h_resultOdd;
      d_result = d_resultOdd;
    }

    size_t chunkElemPairCount;
    if (isLastChunk)
      chunkElemPairCount = elemPairCount - chunk * m_maxActiveElemPairCount;
    else
      chunkElemPairCount = m_maxActiveElemPairCount;

    // Perfrom integration for current element pair chunk on the GPU
    IntegrationTask(taskGroupDevice,
        chunk, isLastChunk, elemPairCount, m_maxActiveElemPairCount,
        d_testElementIndices, d_trialElementIndices,
        m_cudaIntegrators[device], h_result, d_result, resultArraySize,
        integrationTimer);

    // Perfrom assembly for current element pair chunk on the CPU
    AssemblyTask<BasisFunctionType, ResultType, CudaResultType>(taskGroupDevice,
        chunk, isLastChunk, chunkElemPairCount, m_maxActiveElemPairCount,
        m_testDofCount, m_trialDofCount, testElementIndices, trialElementIndices,
        testLocalDofs, trialLocalDofs, testLocalDofWeights, trialLocalDofWeights,
        blockRows, blockCols, h_result, result, mutex, assemblyTimer);
  }

  // Print information about mean integration and assembly time
  if (chunkCount > 2) {

    integrationTimer /= (chunkCount - 2);
    assemblyTimer /= (chunkCount - 2);
    if (integrationTimer > assemblyTimer)
      std::cout << "INFO: Speedup is bound by integration (GPU) with "
          "mean integration time " << integrationTimer << " ms "
          "and mean assembly time " << assemblyTimer << " ms" << std::endl;
    else
      std::cout << "INFO: Speedup is bound by assembly (CPU) with "
          "mean integration time " << integrationTimer << " ms "
          "and mean assembly time " << assemblyTimer << " ms" << std::endl;
  }

  // Wait until last chunk assembly has finished
  taskGroupDevice.wait();

  // Free pageable host memory and memory on the device
  free(h_resultEven);
  cu_verify( cudaFree(d_resultEven) );
  if (chunkCount > 1) {
    free(h_resultOdd);
    cu_verify( cudaFree(d_resultOdd) );
  }
}

} // namespace Fiber
