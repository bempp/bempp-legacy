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

#include "cuda_dense_global_assembler.hpp"

#include "cuda_grid.hpp"
#include "cuda_integrator.hpp"

#include "../assembly/discrete_dense_boundary_operator.hpp"
#include "../assembly/assembly_options.hpp"
#include "../assembly/context.hpp"

#include "../common/types.hpp"
#include "../common/complex_aux.hpp"
#include "../common/global_parameters.hpp"

#include "../fiber/numerical_quadrature.hpp"
#include "../fiber/default_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

#include "../space/space.hpp"

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <chrono>
#include <cuda_profiler_api.h>

namespace Bempp {

// Helper functions and classes
namespace {

template <typename CudaResultType, typename ResultType>
struct AssemblyHelper {
  static ResultType value(
      const CudaResultType *h_regularResult,
      const size_t index, const size_t offset) {
    throw std::invalid_argument(
        "CudaDenseGlobalAssembler::assembleDetachedWeakForm(): "
        "invalid ResultType");
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, float> {
  static float value(
      const CudaResultType *h_regularResult,
      const size_t index, const size_t offset) {
    return static_cast<float>(h_regularResult[index]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, double> {
  static double value(
      const CudaResultType *h_regularResult,
      const size_t index, const size_t offset) {
    return static_cast<double>(h_regularResult[index]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, std::complex<float>> {
  static std::complex<float> value(
      const CudaResultType *h_regularResult,
      const size_t index, const size_t offset) {
    return std::complex<float>(
        h_regularResult[index], h_regularResult[index + offset]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, std::complex<double>> {
  static std::complex<double> value(
      const CudaResultType *h_regularResult,
      const size_t index, const size_t offset) {
    return std::complex<double>(
        h_regularResult[index], h_regularResult[index + offset]);
  }
};

/** Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid() */
template <typename BasisFunctionType>
void gatherGlobalDofs(
    const Space<BasisFunctionType> &space,
    std::vector<std::vector<GlobalDofIndex>> &globalDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) {

  // Get the grid's view so that we can iterate over elements
  const GridView &view = space.gridView();
  const int elementCount = view.entityCount(0);

  // Global DOF indices corresponding to local DOFs on elements
  globalDofs.clear();
  globalDofs.resize(elementCount);
  // Weights of the local DOFs on elements
  localDofWeights.clear();
  localDofWeights.resize(elementCount);

  // Gather global DOF lists
  const Mapper &mapper = view.elementMapper();
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    const int elementIndex = mapper.entityIndex(element);
    space.getGlobalDofs(element, globalDofs[elementIndex],
                        localDofWeights[elementIndex]);
    it->next();
  }
}

// Determine the maximum number of element pairs in a chunk treated on the
// device according to hardware specs
template <typename KernelType, typename CudaCoordinateType>
size_t getMaxActiveElemPairCount(
    const shared_ptr<CudaGrid<CudaCoordinateType>> testGrid,
    const shared_ptr<CudaGrid<CudaCoordinateType>> trialGrid,
    const size_t testIndexCount, const size_t trialIndexCount,
    const unsigned int testDofCount, const unsigned int trialDofCount,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    const CudaOptions cudaOptions, const int deviceId,
    const size_t maxElemPairCountPerDevice, const size_t reserveMem = 1e09) {

  const size_t testElemCount = testGrid->elemCount();
  const size_t trialElemCount = trialGrid->elemCount();

  const unsigned int testDim = testGrid->dim();
  const unsigned int trialDim = trialGrid->dim();

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, deviceId);
  const size_t totalGlobalMem = prop.totalGlobalMem;
  const int warpSize = prop.warpSize;

  size_t gridMem = testElemCount * testGrid->idxCount() * sizeof(int)
        + testGrid->vtxCount() * testDim * sizeof(CudaCoordinateType);
  if (testGrid.get() != trialGrid.get()) {
    gridMem += trialElemCount * trialGrid->idxCount() * sizeof(int)
        + trialGrid->vtxCount() * trialDim * sizeof(CudaCoordinateType);
  }

  const size_t indexMem = (testIndexCount + trialIndexCount) * sizeof(int);

  size_t elemMem = 0;
  if (cudaOptions.isElementDataCachingEnabled()) {
    elemMem +=
        testElemCount * testDim * sizeof(CudaCoordinateType)                      // normals
        + testElemCount * sizeof(CudaCoordinateType);                             // integration elements
        + testElemCount * testDim * testPointCount * sizeof(CudaCoordinateType);  // global quad points
    if (testGrid.get() != trialGrid.get()) {
      elemMem +=
        trialElemCount * trialDim * sizeof(CudaCoordinateType)
        + trialElemCount * sizeof(CudaCoordinateType);
        if (testPointCount != trialPointCount) {
          elemMem +=
              trialElemCount * trialDim * trialPointCount * sizeof(CudaCoordinateType);
        }
    }
  }

  if (totalGlobalMem < (gridMem + indexMem + elemMem + reserveMem))
    throw std::runtime_error(
        "CudaDenseGlobalAssembler::getMaxActiveElemPairCount(): "
        "grid and element data too large");

  size_t maxActiveElemPairCount;
  if (cudaOptions.chunkElemPairCount() == CudaOptions::AUTO) {
    if (cudaOptions.isKernelDataCachingEnabled()) {
      maxActiveElemPairCount =
        (totalGlobalMem - gridMem - indexMem - elemMem - reserveMem) /
        (testPointCount * trialPointCount * sizeof(KernelType));
    } else {
      maxActiveElemPairCount = maxElemPairCountPerDevice;
    }
  } else {
    if (cudaOptions.isKernelDataCachingEnabled()) {
      size_t upperBound  =
          (totalGlobalMem - gridMem - indexMem - elemMem - reserveMem) /
          (testPointCount * trialPointCount * sizeof(KernelType));
      if (cudaOptions.chunkElemPairCount() > upperBound)
        throw std::runtime_error(
            "CudaDenseGlobalAssembler::getMaxActiveElemPairCount(): "
            "chunk size too large");
      maxActiveElemPairCount = cudaOptions.chunkElemPairCount();
    } else {
      maxActiveElemPairCount = cudaOptions.chunkElemPairCount();
    }
  }

  // Let the chunk size be a multiple of the number of test elements
  return std::max((maxActiveElemPairCount / testIndexCount) * testIndexCount,
                  testIndexCount);
  // Let the chunk size be a multiple of the warp size
//  return std::max((maxActiveElemPairCount / warpSize) * warpSize, size_t(warpSize));
}

template <typename BasisFunctionType>
void getParticipatingElements(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    std::vector<int> &testIndices, std::vector<int> &trialIndices,
    std::vector<std::vector<GlobalDofIndex>> &testGlobalDofs,
    std::vector<std::vector<GlobalDofIndex>> &trialGlobalDofs,
    std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
    std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights) {

  // Global DOF indices corresponding to local DOFs on elements
  gatherGlobalDofs(testSpace, testGlobalDofs, testLocalDofWeights);
  if (&testSpace == &trialSpace) {
    trialGlobalDofs = testGlobalDofs;
    trialLocalDofWeights = testLocalDofWeights;
  } else {
    gatherGlobalDofs(trialSpace, trialGlobalDofs, trialLocalDofWeights);
  }
  const size_t testElementCount = testGlobalDofs.size();
  const size_t trialElementCount = trialGlobalDofs.size();

  // Enumerate the test elements that contribute to at least one global DOF
  testIndices.reserve(testElementCount);
  for (size_t testIndex = 0; testIndex < testElementCount; ++testIndex) {
    const int testDofCount = testGlobalDofs[testIndex].size();
    for (int testDof = 0; testDof < testDofCount; ++testDof) {
      int testGlobalDof = testGlobalDofs[testIndex][testDof];
      if (testGlobalDof >= 0) {
        testIndices.push_back(testIndex);
        break;
      }
    }
  }

  // Enumerate the trial elements that contribute to at least one global DOF
  trialIndices.reserve(trialElementCount);
  for (size_t trialIndex = 0; trialIndex < trialElementCount; ++trialIndex) {
    const int trialDofCount = trialGlobalDofs[trialIndex].size();
    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
      if (trialGlobalDof >= 0) {
        trialIndices.push_back(trialIndex);
        break;
      }
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void IntegrationTask(tbb::task_group &taskGroupDevice,
    const size_t chunk, const bool isLastChunk,
    const size_t elemPairOffset, const size_t elemPairCount,
    const size_t maxActiveElemPairCount,
    shared_ptr<Fiber::CudaIntegrator<BasisFunctionType, KernelType, ResultType,
        CudaBasisFunctionType, CudaKernelType, CudaResultType>> &cudaIntegrator,
    CudaResultType *h_regularResult, float &integrationTimer) {

  taskGroupDevice.run_and_wait([
       chunk, isLastChunk, elemPairOffset, elemPairCount, maxActiveElemPairCount,
       &cudaIntegrator, h_regularResult, &integrationTimer] {

    const size_t elemPairIndexBegin = elemPairOffset
        + chunk * maxActiveElemPairCount;

    size_t elemPairIndexEnd;
    if (isLastChunk)
      elemPairIndexEnd = elemPairOffset + elemPairCount;
    else
      elemPairIndexEnd = elemPairOffset + (chunk + 1) * maxActiveElemPairCount;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Evaluate regular integrals over selected element pairs
    cudaIntegrator->integrate(elemPairIndexBegin, elemPairIndexEnd,
        h_regularResult);

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
    const size_t elemPairOffset, const size_t maxActiveElemPairCount,
    const unsigned int testDofCount, const unsigned int trialDofCount,
    const std::vector<int> &testIndices, const std::vector<int> &trialIndices,
    const std::vector<std::vector<GlobalDofIndex>> &testGlobalDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGlobalDofs,
    const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
    const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
    CudaResultType *h_regularResult,
    const Fiber::_2dArray<std::pair<int, Matrix<ResultType>>> &cache,
    Matrix<ResultType> &result, Matrix<tbb::spin_mutex> &mutex,
    float &assemblyTimer) {

  typedef tbb::spin_mutex MutexType;

  taskGroupDevice.run([
      chunk, isLastChunk, chunkElemPairCount,
      elemPairOffset, maxActiveElemPairCount,
      testDofCount, trialDofCount,
      &testIndices, &trialIndices,
      &testGlobalDofs, &trialGlobalDofs,
      &testLocalDofWeights, &trialLocalDofWeights,
      h_regularResult, &cache, &result, &mutex, &assemblyTimer] {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Parallel global assembly
    tbb::parallel_for(size_t(0), size_t(chunkElemPairCount), [
        chunk, isLastChunk, chunkElemPairCount,
        elemPairOffset, maxActiveElemPairCount,
        testDofCount, trialDofCount,
        &testIndices, &trialIndices,
        &testGlobalDofs, &trialGlobalDofs,
        &testLocalDofWeights, &trialLocalDofWeights,
        h_regularResult, &cache, &result, &mutex](size_t chunkElemPair) {

      const size_t offsetResultImag =
          trialDofCount * testDofCount * chunkElemPairCount;

      const size_t offset = elemPairOffset + chunk * maxActiveElemPairCount;

      const size_t testIndexCount = testIndices.size();

      const int trialIndex =
          trialIndices[(offset + chunkElemPair) / testIndexCount];
      const int testIndex =
          testIndices[(offset + chunkElemPair) % testIndexCount];

      // Try to find matrix in singular integral cache
      const Matrix<ResultType> *cachedLocalWeakForm = 0;
      for (size_t n = 0; n < cache.extent(0); ++n) {
        if (cache(n, trialIndex).first == testIndex) {
          cachedLocalWeakForm = &cache(n, trialIndex).second;
          break;
        }
      }

      // Add the integrals to appropriate entries in the operator's matrix
      for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {

        int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
        if (trialGlobalDof < 0)
          continue;
        assert(std::abs(trialLocalDofWeights[trialIndex][trialDof]) > 0.);

        for (int testDof = 0; testDof < testDofCount; ++testDof) {

          int testGlobalDof = testGlobalDofs[testIndex][testDof];
          if (testGlobalDof < 0)
            continue;
          assert(std::abs(testLocalDofWeights[testIndex][testDof]) > 0.);

          MutexType::scoped_lock lock(mutex(testGlobalDof, trialGlobalDof));
//              MutexType::scoped_lock lock(mutex[testGlobalDof]);

          if (cachedLocalWeakForm) { // Matrix found in cache
            result(testGlobalDof, trialGlobalDof) +=
                  conj(testLocalDofWeights[testIndex][testDof]) *
                  trialLocalDofWeights[trialIndex][trialDof] *
                  (*cachedLocalWeakForm)(testDof, trialDof);
          } else {
            const size_t index = trialDof * testDofCount * chunkElemPairCount
                                 + testDof * chunkElemPairCount
                                 + chunkElemPair;
            result(testGlobalDof, trialGlobalDof) +=
                  conj(testLocalDofWeights[testIndex][testDof]) *
                  trialLocalDofWeights[trialIndex][trialDof] *
                  AssemblyHelper<CudaResultType, ResultType>::
                  value(h_regularResult, index, offsetResultImag);
          }
        }
      }
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

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void ParallelDeviceLoop(
    const std::vector<int> &testIndices, const std::vector<int> &trialIndices,
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
    const Fiber::Shapeset<BasisFunctionType> &testShapeset,
    const Fiber::Shapeset<BasisFunctionType> &trialShapeset,
    const Matrix<typename ScalarTraits<BasisFunctionType>::RealType>
    &localTestQuadPoints,
    const Matrix<typename ScalarTraits<BasisFunctionType>::RealType>
    &localTrialQuadPoints,
    const std::vector<typename ScalarTraits<BasisFunctionType>::RealType>
    &testQuadWeights,
    const std::vector<typename ScalarTraits<BasisFunctionType>::RealType>
    &trialQuadWeights,
    const std::vector<std::vector<GlobalDofIndex>> &testGlobalDofs,
    const std::vector<std::vector<GlobalDofIndex>> &trialGlobalDofs,
    const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights,
    const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights,
    const Fiber::_2dArray<std::pair<int, Matrix<ResultType>>> &cache,
    Matrix<ResultType> &result, Matrix<tbb::spin_mutex> &mutex,
    const CudaOptions &cudaOptions) {

  typedef typename ScalarTraits<CudaBasisFunctionType>::RealType CudaCoordinateType;

  // Loop over devices
  tbb::parallel_for (size_t(0), size_t(cudaOptions.devices().size()), [
       &testIndices, &trialIndices,
       &testSpace, &trialSpace,
       &testShapeset, &trialShapeset,
       &localTestQuadPoints, &localTrialQuadPoints,
       &testQuadWeights, &trialQuadWeights,
       &testGlobalDofs, &trialGlobalDofs,
       &testLocalDofWeights, &trialLocalDofWeights,
       &kernel, &cache, &result, &mutex, &cudaOptions](size_t device) {

    const std::vector<int> &deviceIds = cudaOptions.devices();
    const int deviceId = deviceIds[device];
    cu_verify( cudaSetDevice(deviceId) );

    const unsigned int deviceCount = deviceIds.size();

    const unsigned int trialDofCount = trialShapeset.size();
    const unsigned int testDofCount = testShapeset.size();
  //  std::cout << "trialDofCount = " << trialDofCount << ", " << std::flush;
  //  std::cout << "testDofCount = " << testDofCount << std::endl;

    const unsigned int trialPointCount = localTrialQuadPoints.cols();
    const unsigned int testPointCount = localTestQuadPoints.cols();
//    std::cout << "trialPointCount = " << trialPointCount << ", " << std::flush;
//    std::cout << "testPointCount = " << testPointCount << std::endl;

    const size_t trialIndexCount = trialIndices.size();
    const size_t testIndexCount = testIndices.size();
  //  std::cout << "trialIndexCount = " << trialIndexCount
  //            << ", testIndexCount = " << testIndexCount << std::endl;

    const size_t elemPairCountTotal = trialIndexCount * testIndexCount;
  //  std::cout << "elemPairCountTotal = " << elemPairCountTotal << std::endl;

    // Note: This approach assumes identical device types
    const size_t maxElemPairCountPerDevice =
        (elemPairCountTotal - 1) / deviceCount + 1;

    // Get the number of element pairs to be treated on each single device
    size_t elemPairCount = 0;
    if (device == deviceCount-1)
      elemPairCount = elemPairCountTotal - device * maxElemPairCountPerDevice;
    else
      elemPairCount = maxElemPairCountPerDevice;

    // Push raw grid data to the device
    shared_ptr<CudaGrid<CudaCoordinateType>> trialGrid =
        trialSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceId);
    shared_ptr<CudaGrid<CudaCoordinateType>> testGrid =
        testSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceId);

    // Get maximum number of element pairs which can be treated on the device
    // simultaneously
    const size_t maxActiveElemPairCount =
        getMaxActiveElemPairCount<KernelType, CudaCoordinateType>(
            testGrid, trialGrid,
            testIndexCount, trialIndexCount,
            testDofCount, trialDofCount,
            testPointCount, trialPointCount,
            cudaOptions, deviceId, maxElemPairCountPerDevice);

    // Get number of element pair chunks
    const unsigned int chunkCount = static_cast<unsigned int>(
        (elemPairCount-1) / maxActiveElemPairCount + 1);

    // Get element pair offset for the current device
    const size_t elemPairOffset = device * maxElemPairCountPerDevice;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Allocate pinned host memory on the device
    CudaResultType *h_regularResultEven, *h_regularResultOdd;
    size_t size = maxActiveElemPairCount
        * trialDofCount * testDofCount * sizeof(CudaResultType);
    if (ScalarTraits<ResultType>::isComplex) size *= 2;
    cu_verify( cudaHostAlloc((void**)&h_regularResultEven, size,
        cudaHostAllocMapped) );
    cu_verify( cudaHostAlloc((void**)&h_regularResultOdd, size,
        cudaHostAllocMapped) );

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    std::cout << "Time for host vector allocation = "
//              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
//              << " ms" << std::endl;

    // Create CUDA integrator
    shared_ptr<Fiber::CudaIntegrator<
        BasisFunctionType, KernelType, ResultType,
        CudaBasisFunctionType, CudaKernelType, CudaResultType>>
    cudaIntegrator = boost::make_shared<Fiber::CudaIntegrator<
        BasisFunctionType, KernelType, ResultType,
        CudaBasisFunctionType, CudaKernelType, CudaResultType>>(
            localTestQuadPoints, localTrialQuadPoints,
            testQuadWeights, trialQuadWeights,
            testShapeset, trialShapeset,
            testGrid, trialGrid,
            testIndices, trialIndices,
            maxActiveElemPairCount,
            kernel,
            deviceId, cudaOptions);

    tbb::task_group taskGroupDevice;
    float integrationTimer = 0., assemblyTimer = 0.;

    // Loop over chunks of element pairs
    for (size_t chunk = 0; chunk < chunkCount; ++chunk) {

//      std::cout << "chunk = " << chunk << std::endl;

      const bool isLastChunk = (chunk == chunkCount - 1);

      CudaResultType *h_regularResult;
      if (chunk % 2 == 0)
        h_regularResult = h_regularResultEven;
      else
        h_regularResult = h_regularResultOdd;

      size_t chunkElemPairCount;
      if (isLastChunk)
        chunkElemPairCount = elemPairCount - chunk * maxActiveElemPairCount;
      else
        chunkElemPairCount = maxActiveElemPairCount;

      // Perfrom integration for current element pair chunk on the GPU
      IntegrationTask(taskGroupDevice,
          chunk, isLastChunk,
          elemPairOffset, elemPairCount, maxActiveElemPairCount,
          cudaIntegrator, h_regularResult, integrationTimer);

      // Perfrom assembly for current element pair chunk on the CPU
      AssemblyTask(taskGroupDevice,
          chunk, isLastChunk, chunkElemPairCount,
          elemPairOffset, maxActiveElemPairCount,
          testDofCount, trialDofCount,
          testIndices, trialIndices,
          testGlobalDofs, trialGlobalDofs,
          testLocalDofWeights, trialLocalDofWeights,
          h_regularResult, cache, result, mutex, assemblyTimer);
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

    // Free pinned host memory
    cu_verify( cudaFreeHost(h_regularResultEven) );
    cu_verify( cudaFreeHost(h_regularResultOdd) );
  });
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
CudaDenseGlobalAssembler<BasisFunctionType, ResultType>::
    assembleDetachedWeakForm(
        const Space<BasisFunctionType> &testSpace,
        const Space<BasisFunctionType> &trialSpace,
        LocalAssemblerForIntegralOperators &assembler,
        const Context<BasisFunctionType, ResultType> &context) {

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Cast assembler to default assembler
  const Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
      BasisFunctionType, KernelType, ResultType, GeometryFactory>
  &defaultAssembler =
      reinterpret_cast<const Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
          BasisFunctionType, KernelType, ResultType, GeometryFactory> &>(assembler);

  // Get kernel from assembler
  shared_ptr<const Fiber::CollectionOfKernels<KernelType>> kernel =
      defaultAssembler.kernels();

  // Get shapesets from assembler
  shared_ptr<const std::vector<const Shapeset*>>
  testShapesets, trialShapesets;
  defaultAssembler.getShapesets(testShapesets, trialShapesets);
  const Shapeset &testShapeset = *(*testShapesets)[0];
  const Shapeset &trialShapeset = *(*trialShapesets)[0];
  std::cout << "NOTE: Shapesets have to be identical within one space" << std::endl;
  std::cout << "trialDofCount = " << trialShapeset.size() << ", " << std::flush;
  std::cout << "testDofCount = " << testShapeset.size() << std::endl;

  // Get reference to singular integral cache and check if it's active
  const Cache &cache = defaultAssembler.cache();
  bool singularIntegralCachingEnabled =
      context.assemblyOptions().isSingularIntegralCachingEnabled();
  if (!singularIntegralCachingEnabled)
    throw std::invalid_argument(
        "CudaDenseGlobalAssembler::assembleDetachedWeakForm(): "
        "singular integral caching must be enabled");

  // Create a discrete operator represented by a matrix that has to be calculated
  std::unique_ptr<DiscreteDenseBoundaryOperator<ResultType>>
  discreteDenseBoundaryOperator(new DiscreteDenseBoundaryOperator<ResultType>());
  std::cout << "testGlobalDofCount = " << testSpace.globalDofCount()
      << ", trialGlobalDofCount = " << trialSpace.globalDofCount() << std::endl;

  // Create the operator's matrix
  Matrix<ResultType>& result = discreteDenseBoundaryOperator->matrix();
  result.resize(testSpace.globalDofCount(), trialSpace.globalDofCount());
  result.setZero();

  // Calculate the matrix in single or double precision
  if (context.cudaOptions().precision() == "single") {

    typedef typename ScalarTraits<typename ScalarTraits<BasisFunctionType>::RealType>::SingleType CudaBasisFunctionType;
    typedef typename ScalarTraits<typename ScalarTraits<KernelType>::RealType>::SingleType CudaKernelType;
    typedef typename ScalarTraits<typename ScalarTraits<ResultType>::RealType>::SingleType CudaResultType;

    assembleDetachedWeakFormImpl<
    CudaBasisFunctionType, CudaKernelType, CudaResultType>(
        testSpace, trialSpace,
        kernel, testShapeset, trialShapeset, cache,
        context.cudaOptions(), result);

  } else {

    typedef typename ScalarTraits<typename ScalarTraits<BasisFunctionType>::RealType>::DoubleType CudaBasisFunctionType;
    typedef typename ScalarTraits<typename ScalarTraits<KernelType>::RealType>::DoubleType CudaKernelType;
    typedef typename ScalarTraits<typename ScalarTraits<ResultType>::RealType>::DoubleType CudaResultType;

    assembleDetachedWeakFormImpl<
    CudaBasisFunctionType, CudaKernelType, CudaResultType>(
        testSpace, trialSpace,
        kernel, testShapeset, trialShapeset, cache,
        context.cudaOptions(), result);
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for CUDA dense assembly = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  std::ofstream file("cuda_dense_assembly_timer.dat", std::ios::out | std::ios::app);
  file << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

  // Return the discrete operator represented by the matrix that has just been
  // calculated
  return discreteDenseBoundaryOperator;
}

template <typename BasisFunctionType, typename ResultType>
template<typename CudaBasisFunctionType, typename CudaKernelType,
typename CudaResultType>
void CudaDenseGlobalAssembler<BasisFunctionType, ResultType>::
    assembleDetachedWeakFormImpl(
        const Space<BasisFunctionType> &testSpace,
        const Space<BasisFunctionType> &trialSpace,
        const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
        const Shapeset &testShapeset, const Shapeset &trialShapeset,
        const Cache &cache,
        const CudaOptions &cudaOptions,
        Matrix<ResultType> &result) {

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  // Find global DOF indices corresponding to local DOFs on elements and enumerate
  // the test and trial elements that contribute to at least one global DOF
  std::vector<std::vector<GlobalDofIndex>> testGlobalDofs, trialGlobalDofs;
  std::vector<std::vector<BasisFunctionType>> testLocalDofWeights,
    trialLocalDofWeights;
  std::vector<int> testIndices, trialIndices;
  getParticipatingElements(testSpace, trialSpace, testIndices, trialIndices,
      testGlobalDofs, trialGlobalDofs, testLocalDofWeights, trialLocalDofWeights);
  std::cout << "elemPairCountTotal = "
      << trialIndices.size() * testIndices.size() << std::endl;

  // Create a mutex matrix for parallel assembly
  typedef tbb::spin_mutex MutexType;
  Matrix<MutexType> mutex(testSpace.globalDofCount(),
                          trialSpace.globalDofCount());
//  typedef tbb::speculative_spin_mutex MutexType;
//  std::vector<MutexType> mutex(testSpace.globalDofCount());

  // Get numerical quadrature points and weights
  const int trialQuadOrder = cudaOptions.quadOrder();
  const int testQuadOrder = cudaOptions.quadOrder();
  Matrix<CoordinateType> localTrialQuadPoints, localTestQuadPoints;
  std::vector<CoordinateType> trialQuadWeights, testQuadWeights;
  Fiber::fillSingleQuadraturePointsAndWeights(3, trialQuadOrder,
      localTrialQuadPoints, trialQuadWeights);
  Fiber::fillSingleQuadraturePointsAndWeights(3, testQuadOrder,
      localTestQuadPoints, testQuadWeights);
  //  std::cout << "trialQuadOrder = " << trialQuadOrder << ", " << std::flush;
  //  std::cout << "testQuadOrder = " << testQuadOrder << std::endl;

  cudaProfilerStart();

  // Loop over devices to perform matrix computation on multiple GPUs
  ParallelDeviceLoop<BasisFunctionType, KernelType, ResultType,
      CudaBasisFunctionType, CudaKernelType, CudaResultType>(
          testIndices, trialIndices,
          testSpace, trialSpace,
          kernel,
          testShapeset, trialShapeset,
          localTestQuadPoints, localTrialQuadPoints,
          testQuadWeights, trialQuadWeights,
          testGlobalDofs, trialGlobalDofs, testLocalDofWeights, trialLocalDofWeights,
          cache, result, mutex, cudaOptions);

  cudaProfilerStop();
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaDenseGlobalAssembler);

} // namespace Bempp
