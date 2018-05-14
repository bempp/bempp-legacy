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
    CudaResultType *h_result, CudaResultType *d_result,
    const size_t resultArraySize, double &integrationTimer) {

  taskGroupDevice.run_and_wait([
       chunk, isLastChunk, elemPairCount, maxActiveElemPairCount,
       &d_testElementIndices, &d_trialElementIndices,
       &cudaIntegrator, h_result, d_result, resultArraySize, &integrationTimer] {

//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    const size_t elemPairIndexBegin = chunk * maxActiveElemPairCount;

    size_t elemPairIndexEnd;
    if (isLastChunk)
      elemPairIndexEnd = elemPairCount;
    else
      elemPairIndexEnd = (chunk + 1) * maxActiveElemPairCount;

    // Evaluate regular integrals over selected element pairs
    cudaIntegrator->integrate(d_testElementIndices, d_trialElementIndices,
        elemPairIndexBegin, elemPairIndexEnd, d_result);

    // Copy results to host
    cu_verify( cudaMemcpy(h_result, d_result, resultArraySize,
        cudaMemcpyDeviceToHost) );

//    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    integrationTimer += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
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
    double &assemblyTimer) {

//  typedef tbb::spin_mutex MutexType;

  taskGroupDevice.run([
      chunk, isLastChunk, chunkElemPairCount, maxActiveElemPairCount,
      testDofCount, trialDofCount, &testElementIndices, &trialElementIndices,
      &testLocalDofs, &trialLocalDofs, &testLocalDofWeights, &trialLocalDofWeights,
      &blockRows, &blockCols, h_result, &result/*, &mutex*/, &assemblyTimer] {

//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

      // Global assembly
      for (size_t chunkElemPair = 0; chunkElemPair < chunkElemPairCount; ++chunkElemPair) {

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

          result(blockRows[testIndex][testDof],
                 blockCols[trialIndex][trialDof]) +=
              conjugate(testLocalDofWeights[testIndex][testDof]) *
              trialLocalDofWeights[trialIndex][trialDof] *
              AssemblyHelper<CudaResultType, ResultType>::
              value(h_result, index, offsetResultImag);
        }
      }
    };
//    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    assemblyTimer += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

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
        const Bempp::Context<BasisFunctionType, ResultType> &context)
        : m_testDofCount(testShapeset.size()),
          m_trialDofCount(trialShapeset.size()) {

  if (m_trialDofCount != 3 || m_testDofCount != 3)
    throw std::invalid_argument(
        "CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces::CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces(): "
        "Only three local dofs per element supported on the device");

  // Get numerical quadrature points and weights
  int quadOrder[3];
  quadOrder[0] = context.globalParameterList().template get<int>("options.quadrature.near.doubleOrder");
  quadOrder[1] = context.globalParameterList().template get<int>("options.quadrature.medium.doubleOrder");
  quadOrder[2] = context.globalParameterList().template get<int>("options.quadrature.far.doubleOrder");
  if (quadOrder[0] != 4)
    throw std::invalid_argument(
        "CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces::CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces(): "
        "Maximum of six quadrature points per element supported on the device");
  std::vector<      Matrix<CoordinateType> > localTrialQuadPoints(3), localTestQuadPoints(3);
  std::vector< std::vector<CoordinateType> >     trialQuadWeights(3),     testQuadWeights(3);
  for (int q = 0; q < 3; ++q) {
    fillSingleQuadraturePointsAndWeights(3, quadOrder[q],
        localTrialQuadPoints[q], trialQuadWeights[q]);
    fillSingleQuadraturePointsAndWeights(3, quadOrder[q],
        localTestQuadPoints[q], testQuadWeights[q]);
  }

  m_deviceIds = context.cudaOptions().devices();
  const unsigned int deviceCount = m_deviceIds.size();

  m_cudaIntegrators.resize(deviceCount);
  m_buffer.resize(deviceCount);

  // Get maximum number of element pairs which can be treated on the device
  // simultaneously
  m_maxActiveElemPairCount = context.cudaOptions().chunkElemPairCount();

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
        m_maxActiveElemPairCount, kernel, deviceId, context.cudaOptions());

    m_buffer[device] = tbb::enumerable_thread_specific<
        ::CudaHostDeviceBuffer<ResultType, CudaResultType> >(deviceId, m_maxActiveElemPairCount,
            m_trialDofCount, m_testDofCount);
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
        Matrix<ResultType> &result, const unsigned int device/*,
        tbb::atomic<std::chrono::steady_clock::duration>& allocationTimer,
        tbb::atomic<std::chrono::steady_clock::duration>& integrationTimer,
        tbb::atomic<std::chrono::steady_clock::duration>& assemblyTimer*/) {

  double localTotalTimer = 0., localAllocationTimer = 0., localIntegrationTimer = 0., localAssemblyTimer = 0.;
//  std::chrono::steady_clock::time_point totalBegin = std::chrono::steady_clock::now();

  const int deviceId = m_deviceIds[device];
  cu_verify( cudaSetDevice(deviceId) );

  // Get the number of element pairs to be treated on the specified device
  const size_t elemPairCount = trialElementIndices.size() * testElementIndices.size();

  // Get number of element pair chunks
  const unsigned int chunkCount = static_cast<unsigned int>(
      (elemPairCount-1) / m_maxActiveElemPairCount + 1);

//  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  CudaResultType *h_resultEven, *h_resultOdd, *d_resultEven, *d_resultOdd;
  size_t resultArraySize;
  if (chunkCount == 1) {

    resultArraySize = elemPairCount
        * m_trialDofCount * m_testDofCount * sizeof(CudaResultType);
    if (Bempp::ScalarTraits<ResultType>::isComplex) resultArraySize *= 2;

    // Allocate pageable host memory and memory on the device
    h_resultEven = (CudaResultType*)malloc(resultArraySize);
    cu_verify( cudaMalloc((void**)&(d_resultEven), resultArraySize) );

  } else {

    resultArraySize = m_maxActiveElemPairCount
        * m_trialDofCount * m_testDofCount * sizeof(CudaResultType);
    if (Bempp::ScalarTraits<ResultType>::isComplex) resultArraySize *= 2;

    // Allocate pageable host memory and memory on the device
    h_resultEven = (CudaResultType*)malloc(resultArraySize);
    h_resultOdd = (CudaResultType*)malloc(resultArraySize);
    cu_verify( cudaMalloc((void**)&(d_resultEven), resultArraySize) );
    cu_verify( cudaMalloc((void**)&(d_resultOdd), resultArraySize) );
  }

//  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//  localAllocationTimer = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
//  std::cout << "Time for memory allocation = " << localAllocationTimer << " us" << std::endl;

  // Copy element indices to device memory
  thrust::device_vector<int> d_testElementIndices, d_trialElementIndices;
  m_cudaIntegrators[device]->pushElementIndicesToDevice(
      testElementIndices, trialElementIndices,
      d_testElementIndices, d_trialElementIndices);

  tbb::task_group taskGroupDevice;

  // Loop over chunks of element pairs
  for (size_t chunk = 0; chunk < chunkCount; ++chunk) {

    const bool isLastChunk = (chunk == chunkCount - 1);

    CudaResultType *h_result, *d_result;
    if (chunk % 2 == 0) {
      h_result = h_resultEven;
      d_result = d_resultEven;
//      h_result = m_buffer[device].local().h_resultEven;
//      d_result = m_buffer[device].local().d_resultEven;
    } else {
      h_result = h_resultOdd;
      d_result = d_resultOdd;
//      h_result = m_buffer[device].local().h_resultOdd;
//      d_result = m_buffer[device].local().d_resultOdd;
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
        localIntegrationTimer);

    // Perfrom assembly for current element pair chunk on the CPU
    AssemblyTask<BasisFunctionType, ResultType, CudaResultType>(taskGroupDevice,
        chunk, isLastChunk, chunkElemPairCount, m_maxActiveElemPairCount,
        m_testDofCount, m_trialDofCount, testElementIndices, trialElementIndices,
        testLocalDofs, trialLocalDofs, testLocalDofWeights, trialLocalDofWeights,
        blockRows, blockCols, h_result, result, localAssemblyTimer);
  }

  // Wait until last chunk assembly has finished
  taskGroupDevice.wait();

//  std::cout << "Time for integration = " << localIntegrationTimer << " us" << std::endl;
//  std::cout << "Time for assembly = "    << localAssemblyTimer    << " us" << std::endl;

  // Free pageable host memory and memory on the device
  free(h_resultEven);
  cu_verify( cudaFree(d_resultEven) );
  if (chunkCount > 1) {
    free(h_resultOdd);
    cu_verify( cudaFree(d_resultOdd) );
  }

//  std::chrono::steady_clock::time_point totalEnd = std::chrono::steady_clock::now();
//  localTotalTimer = std::chrono::duration_cast<std::chrono::microseconds>(totalEnd - totalBegin).count();
//  std::cout << "Time for total block computation = " << localTotalTimer << " us" << std::endl;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaDefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    evaluateLocalWeakForms(
        const std::vector<std::vector<Bempp::LocalDof>> & testRawLocalDofs,
        const std::vector<std::vector<Bempp::LocalDof>> &trialRawLocalDofs,
        Matrix<ResultType> &result, const unsigned int device) {

  typedef std::chrono::steady_clock clock;
  double totalTimer = 0., allocTimer = 0., h2dTimer = 0., intTimer = 0., d2hTimer = 0.;
  clock::time_point totalBegin = clock::now();

  const int deviceId = m_deviceIds[device];
  cu_verify( cudaSetDevice(deviceId) );

  // Allocate device memory
  CudaResultType     *h_result              , *d_result;
  Bempp::EntityIndex *h_testElemIndices     , *d_testElemIndices;
  char               *h_testLocalDofIndices , *d_testLocalDofIndices;
  Bempp::EntityIndex *h_trialElemIndices    , *d_trialElemIndices;
  char               *h_trialLocalDofIndices, *d_trialLocalDofIndices;

  const size_t maxElemsPerDof = 10; // TODO Check that!
  size_t resultArraySize;
  resultArraySize = testRawLocalDofs.size() * trialRawLocalDofs.size() * sizeof(CudaResultType);
  if (Bempp::ScalarTraits<ResultType>::isComplex) resultArraySize *= 2;

  typedef typename std::conditional<Bempp::ScalarTraits<ResultType>::isComplex,
                           typename Bempp::ScalarTraits<CudaResultType>::ComplexType,
                           CudaResultType>::type CudaComplexResultType;

  clock::time_point allocBegin = clock::now();
  Matrix<CudaComplexResultType> cudaComplexResult(result.rows(), result.cols());
  h_result = reinterpret_cast<CudaResultType*>(cudaComplexResult.data());
  cu_verify( cudaMalloc((void**)&(d_result), resultArraySize) );
  cu_verify( cudaMalloc((void**)&(d_testElemIndices     ), testRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex)) );
  cu_verify( cudaMalloc((void**)&(d_testLocalDofIndices ), testRawLocalDofs.size() * maxElemsPerDof * sizeof(char              )) );
  h_testElemIndices     = (Bempp::EntityIndex*)malloc(     testRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex))  ;
  h_testLocalDofIndices = (char              *)malloc(     testRawLocalDofs.size() * maxElemsPerDof * sizeof(char              ))  ;
  cu_verify( cudaMalloc((void**)&(d_trialElemIndices    ),trialRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex)) );
  cu_verify( cudaMalloc((void**)&(d_trialLocalDofIndices),trialRawLocalDofs.size() * maxElemsPerDof * sizeof(char              )) );
  h_trialElemIndices     = (Bempp::EntityIndex*)malloc(   trialRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex))  ;
  h_trialLocalDofIndices = (char              *)malloc(   trialRawLocalDofs.size() * maxElemsPerDof * sizeof(char              ))  ;
  clock::time_point allocEnd = clock::now();
  allocTimer = std::chrono::duration_cast<std::chrono::microseconds>(allocEnd - allocBegin).count();

  // Gather dof information and copy to device
  memset(h_testElemIndices, Bempp::EntityIndex(-1), testRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex));
  for (int i = 0; i < testRawLocalDofs.size(); ++i) {
    assert(testRawLocalDofs[i].size() < maxElemsPerDof);
    for (int j = 0; j < testRawLocalDofs[i].size(); ++j) {
      h_testElemIndices    [i * maxElemsPerDof + j] = testRawLocalDofs[i][j].entityIndex;
      h_testLocalDofIndices[i * maxElemsPerDof + j] = testRawLocalDofs[i][j].dofIndex;
    }
  }
  memset(h_trialElemIndices, Bempp::EntityIndex(-1), trialRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex));
  for (int i = 0; i < trialRawLocalDofs.size(); ++i) {
    assert(trialRawLocalDofs[i].size() < maxElemsPerDof);
    for (int j = 0; j < trialRawLocalDofs[i].size(); ++j) {
      h_trialElemIndices    [i * maxElemsPerDof + j] = trialRawLocalDofs[i][j].entityIndex;
      h_trialLocalDofIndices[i * maxElemsPerDof + j] = trialRawLocalDofs[i][j].dofIndex;
    }
  }

  clock::time_point h2dBegin = clock::now();
  cu_verify( cudaMemcpyAsync((void*)d_testElemIndices     , (const void*)h_testElemIndices     , testRawLocalDofs.size()  * maxElemsPerDof * sizeof(Bempp::EntityIndex), cudaMemcpyHostToDevice, cudaStreamPerThread) );
  cu_verify( cudaMemcpyAsync((void*)d_testLocalDofIndices , (const void*)h_testLocalDofIndices , testRawLocalDofs.size()  * maxElemsPerDof * sizeof(char              ), cudaMemcpyHostToDevice, cudaStreamPerThread) );
  cu_verify( cudaMemcpyAsync((void*)d_trialElemIndices    , (const void*)h_trialElemIndices    , trialRawLocalDofs.size() * maxElemsPerDof * sizeof(Bempp::EntityIndex), cudaMemcpyHostToDevice, cudaStreamPerThread) );
  cu_verify( cudaMemcpyAsync((void*)d_trialLocalDofIndices, (const void*)h_trialLocalDofIndices, trialRawLocalDofs.size() * maxElemsPerDof * sizeof(char              ), cudaMemcpyHostToDevice, cudaStreamPerThread) );
  clock::time_point h2dEnd = clock::now();
  h2dTimer = std::chrono::duration_cast<std::chrono::microseconds>(h2dEnd - h2dBegin).count();

  // Evaluate coefficients from regular integrals for requested dofs
  clock::time_point intBegin = clock::now();
  m_cudaIntegrators[device]->integrate(
      d_testElemIndices,  d_testLocalDofIndices,  testRawLocalDofs.size(),
      d_trialElemIndices, d_trialLocalDofIndices, trialRawLocalDofs.size(),
      d_result);
  clock::time_point intEnd = clock::now();
  intTimer = std::chrono::duration_cast<std::chrono::microseconds>(intEnd - intBegin).count();

  // Copy results to host memory and cast to host data type
  clock::time_point d2hBegin = clock::now();
  cu_verify( cudaMemcpyAsync(h_result, d_result, resultArraySize, cudaMemcpyDeviceToHost, cudaStreamPerThread) );
  result = cudaComplexResult.template cast<ResultType>();
  clock::time_point d2hEnd = clock::now();
  d2hTimer = std::chrono::duration_cast<std::chrono::microseconds>(d2hEnd - d2hBegin).count();

  // Free memory on host and device
  free(h_testElemIndices     ); cu_verify( cudaFree(d_testElemIndices     ) );
  free(h_testLocalDofIndices ); cu_verify( cudaFree(d_testLocalDofIndices ) );
  free(h_trialElemIndices    ); cu_verify( cudaFree(d_trialElemIndices    ) );
  free(h_trialLocalDofIndices); cu_verify( cudaFree(d_trialLocalDofIndices) );

  clock::time_point totalEnd = clock::now();
  totalTimer = std::chrono::duration_cast<std::chrono::microseconds>(totalEnd - totalBegin).count();
//  printf("N_dofs = %d, t_total = %8.1f us, t_alloc = %8.1f us, t_h2d = %8.1f, t_int = %8.1f us, t_d2h = %8.1f.\n",
//      testRawLocalDofs.size() * trialRawLocalDofs.size(),
//      totalTimer, allocTimer, h2dTimer, intTimer, d2hTimer);
}

} // namespace Fiber
