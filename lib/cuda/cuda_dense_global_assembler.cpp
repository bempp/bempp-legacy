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

#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/numerical_quadrature.hpp"
#include "../fiber/quadrature_descriptor_selector_for_integral_operators.hpp"
#include "../fiber/double_quadrature_rule_family.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/default_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

#include "../space/space.hpp"
#include "../grid/grid.hpp"

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <algorithm>
#include <numeric>
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
        h_regularResult[index], h_regularResult[index+offset]);
  }
};

template <typename CudaResultType>
struct AssemblyHelper<CudaResultType, std::complex<double>> {
  static std::complex<double> value(
      const CudaResultType *h_regularResult,
      const size_t index, const size_t offset) {
    return std::complex<double>(
        h_regularResult[index], h_regularResult[index+offset]);
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

  // Let the chunk size be a multiple of the warp size
//  return (maxActiveElemPairCount / warpSize) * warpSize;

  // Let the chunk size be a multiple of the number of test elements
  return (maxActiveElemPairCount / testIndexCount) * testIndexCount;
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

  // Create the operator's matrix
  Matrix<ResultType>& result = discreteDenseBoundaryOperator->matrix();
  result.resize(testSpace.globalDofCount(), trialSpace.globalDofCount());
  result.setZero();
    std::cout << "testGlobalDofCount = " << testSpace.globalDofCount()
        << ", trialGlobalDofCount = " << trialSpace.globalDofCount() << std::endl;

  // Calculate the matrix
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
  typedef typename ScalarTraits<CudaBasisFunctionType>::RealType CudaCoordinateType;

  // Global DOF indices corresponding to local DOFs on elements
  std::vector<std::vector<GlobalDofIndex>> testGlobalDofs, trialGlobalDofs;
  std::vector<std::vector<BasisFunctionType>> testLocalDofWeights,
    trialLocalDofWeights;
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
  std::vector<int> testIndices;
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
  std::vector<int> trialIndices;
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

  // TODO: Find a smaller mutex data type
  typedef tbb::spin_mutex MutexType;
  Matrix<MutexType> mutex(testSpace.globalDofCount(),
                          trialSpace.globalDofCount());
//  typedef tbb::speculative_spin_mutex MutexType;
//  std::vector<MutexType> mutex(testSpace.globalDofCount());

  const std::vector<int> &deviceIds = cudaOptions.devices();
  const unsigned int deviceCount = deviceIds.size();

  const size_t trialIndexCount = trialIndices.size();
  const size_t testIndexCount = testIndices.size();
//  std::cout << "trialIndexCount = " << trialIndexCount
//            << ", testIndexCount = " << testIndexCount << std::endl;

  const size_t elemPairCountTotal = trialIndexCount * testIndexCount;
//  std::cout << "elemPairCountTotal = " << elemPairCountTotal << std::endl;

  // Note: This approach assumes identical device types
  const size_t maxElemPairCountPerDevice = (elemPairCountTotal-1) / deviceCount + 1;

  // Get number of element pairs for all active devices
  std::vector<size_t> elemPairCounts(deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
    if (device == deviceCount-1) {
      elemPairCounts[device] = elemPairCountTotal - device * maxElemPairCountPerDevice;
    } else {
      elemPairCounts[device] = maxElemPairCountPerDevice;
    }
  }
//  std::cout << "elemPairCounts = " << std::endl;
//  for (int i = 0; i < elemPairCounts.size(); ++i) {
//    std::cout << elemPairCounts[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  cudaProfilerStart();

  // TODO: Push raw grid data to all active devices (maybe in parallel)
  std::vector<shared_ptr<CudaGrid<CudaCoordinateType>>> trialGrids(deviceCount);
  std::vector<shared_ptr<CudaGrid<CudaCoordinateType>>> testGrids(deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
    trialGrids[device] =
        trialSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceIds[device]);
    testGrids[device] =
        testSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceIds[device]);
  }

  const unsigned int trialDofCount = trialShapeset.size();
  const unsigned int testDofCount = testShapeset.size();
//  std::cout << "trialDofCount = " << trialDofCount << ", " << std::flush;
//  std::cout << "testDofCount = " << testDofCount << std::endl;

  // Get the quadrature order from CUDA options
  const int trialQuadOrder = cudaOptions.quadOrder();
  const int testQuadOrder = cudaOptions.quadOrder();
//  std::cout << "trialQuadOrder = " << trialQuadOrder << ", " << std::flush;
//  std::cout << "testQuadOrder = " << testQuadOrder << std::endl;
  Matrix<CoordinateType> localTrialQuadPoints, localTestQuadPoints;
  std::vector<CoordinateType> trialQuadWeights, testQuadWeights;
  Fiber::fillSingleQuadraturePointsAndWeights(
      trialGrids[0]->idxCount(), trialQuadOrder,
      localTrialQuadPoints, trialQuadWeights);
  Fiber::fillSingleQuadraturePointsAndWeights(
      testGrids[0]->idxCount(), testQuadOrder,
      localTestQuadPoints, testQuadWeights);
  const unsigned int trialPointCount = localTrialQuadPoints.cols();
  const unsigned int testPointCount = localTestQuadPoints.cols();

  // Get maximum number of active element pairs for all active devices
  std::vector<size_t> maxActiveElemPairCounts(deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
    maxActiveElemPairCounts[device] =
        getMaxActiveElemPairCount<KernelType, CudaCoordinateType>(
            testGrids[device], trialGrids[device],
            testIndexCount, trialIndexCount,
            testDofCount, trialDofCount,
            testPointCount, trialPointCount,
            cudaOptions, deviceIds[device], maxElemPairCountPerDevice);
  }
//  std::cout << "maxActiveElemPairCounts = " << std::endl;
//  for (int i = 0; i < maxActiveElemPairCounts.size(); ++i) {
//    std::cout << maxActiveElemPairCounts[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  // Get number of chunks for all active devices
  std::vector<unsigned int> chunkCounts(deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
    chunkCounts[device] = static_cast<unsigned int>(
            (elemPairCounts[device]-1) / maxActiveElemPairCounts[device] + 1);
  }
//  std::cout << "chunkCounts = " << std::endl;
//  for (int i = 0; i < chunkCounts.size(); ++i) {
//    std::cout << chunkCounts[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  // Allocate pinned host memory on all active devices
  std::vector<CudaResultType*> h_regularResultEven(deviceCount);
  std::vector<CudaResultType*> h_regularResultOdd(deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
    size_t size = maxActiveElemPairCounts[device]
        * trialDofCount * testDofCount * sizeof(CudaResultType);
    if (ScalarTraits<ResultType>::isComplex) size *= 2;
    cudaSetDevice(deviceIds[device]);
    cu_verify( cudaHostAlloc((void**)&h_regularResultEven[device], size,
        cudaHostAllocMapped) );
    cu_verify( cudaHostAlloc((void**)&h_regularResultOdd[device], size,
        cudaHostAllocMapped) );
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//  std::cout << "Time for host vector allocation = "
//            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
//            << " ms" << std::endl;

  // Create CUDA integrators on all active devices
  std::vector<shared_ptr<Fiber::CudaIntegrator<
      BasisFunctionType, KernelType, ResultType,
      CudaBasisFunctionType, CudaKernelType, CudaResultType>>>
  cudaIntegrators(deviceCount);
  for (int device = 0; device < deviceCount; ++device) {
    cudaIntegrators[device] = boost::make_shared<Fiber::CudaIntegrator<
        BasisFunctionType, KernelType, ResultType,
        CudaBasisFunctionType, CudaKernelType, CudaResultType>>(
            localTestQuadPoints, localTrialQuadPoints,
            testQuadWeights, trialQuadWeights,
            testShapeset, trialShapeset,
            testGrids[device], trialGrids[device],
            testIndices, trialIndices,
            kernel,
            deviceIds[device], cudaOptions);
  }

  // Loop over devices
  tbb::parallel_for (size_t(0), size_t(deviceCount),[
       &chunkCounts, &elemPairCounts, &maxActiveElemPairCounts, &deviceIds,
       testDofCount, trialDofCount,
       &testIndices, &trialIndices,
       &cudaIntegrators, &h_regularResultEven, &h_regularResultOdd,
       &testGlobalDofs, &trialGlobalDofs,
       &testLocalDofWeights, &trialLocalDofWeights,
       cache, &result, &mutex](size_t device) {

    double integration_timer = 0., assembly_timer = 0.;

    const unsigned int chunkCount = chunkCounts[device];
    size_t elemPairOffset = 0;
    for (int i = 0; i < device; ++i) {
      elemPairOffset += elemPairCounts[i];
    }

    tbb::task_group taskGroupDevice;

    // Loop over element pair chunks
    for (int chunk = 0; chunk < chunkCount; ++chunk) {

//      std::cout << "chunk = " << chunk << std::endl;

      CudaResultType *h_regularResult;
      if (chunk % 2 == 0)
        h_regularResult = h_regularResultEven[device];
      else
        h_regularResult = h_regularResultOdd[device];

      size_t chunkElemPairCount;
      if (chunk == chunkCount-1) {
        chunkElemPairCount = elemPairCounts[device]
            - chunk * maxActiveElemPairCounts[device];
      } else {
        chunkElemPairCount = maxActiveElemPairCounts[device];
      }

      taskGroupDevice.run_and_wait([device, chunk, chunkCount,
           elemPairOffset, &elemPairCounts, &maxActiveElemPairCounts,
           chunkElemPairCount, testDofCount, trialDofCount,
           &cudaIntegrators, h_regularResult, &integration_timer]{

        const size_t elemPairIndexBegin = elemPairOffset
            + chunk * maxActiveElemPairCounts[device];
        size_t elemPairIndexEnd;
        if (chunk == chunkCount-1) {
          elemPairIndexEnd = elemPairOffset + elemPairCounts[device];
        } else {
          elemPairIndexEnd = elemPairOffset
              + (chunk+1) * maxActiveElemPairCounts[device];
        }

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        // Evaluate regular integrals over selected element pairs
        cudaIntegrators[device]->integrate(
            elemPairIndexBegin, elemPairIndexEnd,
            h_regularResult);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        integration_timer += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        std::cout << "Time for CudaIntegrator::integrate() = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                  << " ms" << std::endl;
      });

      taskGroupDevice.run([device, chunk, chunkCount,
          elemPairOffset, &maxActiveElemPairCounts,
          chunkElemPairCount, testDofCount, trialDofCount,
          &testIndices, &trialIndices,
          &testGlobalDofs, &trialGlobalDofs,
          &testLocalDofWeights, &trialLocalDofWeights,
          h_regularResult, cache, &result, &mutex, &assembly_timer]{

        const size_t offset = elemPairOffset
            + chunk * maxActiveElemPairCounts[device];
        const size_t testIndexCount = testIndices.size();

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        // Global assembly
        tbb::parallel_for(size_t(0), size_t(chunkElemPairCount), [
            device, chunkElemPairCount, offset, testIndexCount,
            trialDofCount, testDofCount,
            &testIndices, &trialIndices,
            &testGlobalDofs, &trialGlobalDofs,
            &testLocalDofWeights, &trialLocalDofWeights,
            h_regularResult, cache, &result, &mutex](size_t chunkElemPair) {

          const size_t offsetResultImag =
              trialDofCount * testDofCount * chunkElemPairCount;

          const int trialIndex =
              trialIndices[(offset + chunkElemPair) / testIndexCount];
          const int testIndex =
              testIndices[(offset + chunkElemPair) % testIndexCount];

          // Try to find matrix in singular integral cache
          const Matrix<ResultType> *cachedLocalWeakForm = 0;
          for (size_t n = 0; n < cache.extent(0); ++n)
            if (cache(n, trialIndex).first == testIndex) {
              cachedLocalWeakForm = &cache(n, trialIndex).second;
              break;
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
        assembly_timer += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        std::cout << "Time for regular global result assembly = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                  << " ms" << std::endl;
      });
    }
    integration_timer /= chunkCount;
    assembly_timer /= chunkCount;
    if (integration_timer > assembly_timer)
      std::cout << "INFO: Speedup is bound by integration (GPU) with "
          "mean integration time " << integration_timer << " ms "
          "and mean assembly time " << assembly_timer << " ms" << std::endl;
    else
      std::cout << "INFO: Speedup is bound by assembly (CPU) with "
          "mean integration time " << integration_timer << " ms "
          "and mean assembly time " << assembly_timer << " ms" << std::endl;

    taskGroupDevice.wait();
  });
  cudaProfilerStop();

  if (result.rows() < 25 && result.cols() < 25) {
    std::cout << "result (cudadense) = " << std::endl;
    std::cout << result << std::endl;
  }

  // Free pinned host memory
  for (int device = 0; device < deviceCount; ++device) {
    cudaSetDevice(deviceIds[device]);
    cu_verify( cudaFreeHost(h_regularResultEven[device]) );
    cu_verify( cudaFreeHost(h_regularResultOdd[device]) );
  }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaDenseGlobalAssembler);

} // namespace Bempp
