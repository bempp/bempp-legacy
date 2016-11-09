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
#include "cuda_options.hpp"

#include "../assembly/discrete_dense_boundary_operator.hpp"
#include "../assembly/assembly_options.hpp"
#include "../assembly/context.hpp"

#include "../common/types.hpp"
#include "../common/complex_aux.hpp"
#include "../common/global_parameters.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/element_pair_topology.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/numerical_quadrature.hpp"
#include "../fiber/quadrature_descriptor_selector_for_integral_operators.hpp"
#include "../fiber/double_quadrature_rule_family.hpp"
#include "../fiber/serial_blas_region.hpp"

#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

#include "../space/space.hpp"
#include "../grid/grid.hpp"

#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <tbb/queuing_mutex.h>
#include <tbb/concurrent_vector.h>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cuda_profiler_api.h>

namespace Bempp {

// Helper functions and classes
namespace {

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
template <typename ResultType>
unsigned int getMaxActiveElemPairCount(
    shared_ptr<CudaGrid<typename ScalarTraits<ResultType>::RealType>> testGrid,
    shared_ptr<CudaGrid<typename ScalarTraits<ResultType>::RealType>> trialGrid,
    const unsigned int testDofCount, const unsigned int trialDofCount,
    const unsigned int testPointCount, const unsigned int trialPointCount,
    bool cacheElemData, size_t reserveMem = 1e09) {

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  const unsigned int testElemCount = testGrid->elemCount();
  const unsigned int trialElemCount = trialGrid->elemCount();

  const unsigned int testDim = testGrid->dim();
  const unsigned int trialDim = trialGrid->dim();

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  size_t totalGlobalMem = prop.totalGlobalMem;
  const int warpSize = prop.warpSize;

  size_t gridMem = testElemCount * testGrid->idxCount() * sizeof(int)
        + testGrid->vtxCount() * testDim * sizeof(CoordinateType);
  if (testGrid.get() != trialGrid.get()) {
    gridMem += trialElemCount * trialGrid->idxCount() * sizeof(int)
        + trialGrid->vtxCount() * trialDim * sizeof(CoordinateType);
  }

  size_t elemMem = 0;
  if (cacheElemData == true) {
    elemMem +=
        testElemCount * testDim * sizeof(CoordinateType) // normals
        + testElemCount * sizeof(CoordinateType);        // integration elements
        + testElemCount * testDim * testPointCount;      // global quad points
    if (testGrid.get() != trialGrid.get() || testPointCount != trialPointCount) {
      elemMem +=
      trialElemCount * trialDim * sizeof(CoordinateType)
      + trialElemCount * sizeof(CoordinateType);
      + trialElemCount * trialDim * trialPointCount;
    }
  }

  int maxActiveElemPairCount =
      (totalGlobalMem - gridMem - elemMem - reserveMem) / (2 * sizeof(int)
      + testDofCount * trialDofCount * sizeof(ResultType));

  // Let the chunk size be a multiple of the warp size
  return (maxActiveElemPairCount / warpSize) * warpSize;
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
CudaDenseGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    LocalAssemblerForIntegralOperators &assembler,
    const Context<BasisFunctionType, ResultType> &context) {

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

  const int testElementCount = testGlobalDofs.size();
  const int trialElementCount = trialGlobalDofs.size();

  // Enumerate the test elements that contribute to at least one global DOF
  std::vector<int> testIndices;
  testIndices.reserve(testElementCount);

  for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
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

  for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex) {
    const int trialDofCount = trialGlobalDofs[trialIndex].size();
    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
      int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
      if (trialGlobalDof >= 0) {
        trialIndices.push_back(trialIndex);
        break;
      }
    }
  }

//  std::cout << "testIndices = " << std::endl;
//  for (int i = 0; i < testIndices.size(); ++i) {
//    std::cout << testIndices[i] << " " << std::flush;
//  }
//  std::cout << std::endl;
//
//  std::cout << "trialIndices = " << std::endl;
//  for (int i = 0; i < trialIndices.size(); ++i) {
//    std::cout << trialIndices[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  // Create the operator's matrix
  Matrix<ResultType> result(testSpace.globalDofCount(),
                            trialSpace.globalDofCount());
  result.setZero();

  tbb::task_group taskGroupGlobal;
  // TODO: Most suitable mutex type?
  tbb::spin_mutex mutex;

  // TODO: How to determine the correct kernel type?
  typedef CoordinateType KernelType;         // Real kernel
//        typedef ResultType KernelType;             // Complex kernel
  shared_ptr<const Fiber::CollectionOfKernels<KernelType>> kernels;
//        assembler.getKernels(kernels);

  // Define CUDA operation settings
  Fiber::CudaOptions cudaOptions;
  cudaOptions.setStreamCount(4);
  cudaOptions.setDevices({0,1});

  const std::vector<int> &devices = cudaOptions.devices();
  const unsigned int deviceCount = devices.size();

  cudaProfilerStart();

  std::vector<shared_ptr<CudaGrid<CoordinateType>>> testGrids(deviceCount);
  std::vector<shared_ptr<CudaGrid<CoordinateType>>> trialGrids(deviceCount);
  // TODO: Push raw grid data to devices
  for (int device = 0; device < deviceCount; ++device) {
    const int deviceId = devices[device];
    testGrids[device] =
        testSpace.grid()->template pushToDevice<CoordinateType>(deviceId);
    trialGrids[device] =
        trialSpace.grid()->template pushToDevice<CoordinateType>(deviceId);
  }
  shared_ptr<CudaGrid<CoordinateType>> testGrid =
      testSpace.grid()->template pushToDevice<CoordinateType>(0);
  shared_ptr<CudaGrid<CoordinateType>> trialGrid =
      trialSpace.grid()->template pushToDevice<CoordinateType>(0);

  // TODO: Is it really necessary to specify the element pair indices?
  const unsigned int testIndexCount = testIndices.size();
  const unsigned int trialIndexCount = trialIndices.size();
  std::vector<int> testElemPairIndices(testIndexCount * trialIndexCount);
  std::vector<int> trialElemPairIndices(testIndexCount * trialIndexCount);
  for (int i = 0; i < testIndexCount; ++i) {
    const int testIndex = testIndices[i];
    for (int j = 0; j < trialIndexCount; ++j) {
      const int trialIndex = trialIndices[j];
      testElemPairIndices[i * trialIndexCount + j] = testIndex;
      trialElemPairIndices[i * trialIndexCount + j] = trialIndex;
    }
  }

  const unsigned int elemPairCount = testElemPairIndices.size();
  if (elemPairCount != trialElemPairIndices.size())
    throw std::invalid_argument(
        "CudaDenseGlobalAssembler::assembleDetachedWeakForm(): "
        "numbers of test and trial elements do not match");
  std::cout << "elemPairCount = " << elemPairCount << std::endl;

  // TODO: Get it from assembler
  std::vector<const Shapeset*> testShapesets, trialShapesets;
  getAllShapesets(testSpace, testShapesets);
  if (&testSpace == &trialSpace) {
    trialShapesets = testShapesets;
  } else {
    getAllShapesets(trialSpace, trialShapesets);
  }

  // Note: Shapesets have to be identical within one space!
  const Shapeset &testShapeset = *testShapesets[0];
  const Shapeset &trialShapeset = *trialShapesets[0];

  const unsigned int testDofCount = testShapeset.size();
  const unsigned int trialDofCount = trialShapeset.size();
  std::cout << "testDofCount = " << testDofCount << ", " << std::flush;
  std::cout << "trialDofCount = " << trialDofCount << std::endl;

  // TODO: Get the quadrature order from CUDA options
  const auto parameterList = context.globalParameterList();
  auto nearDoubleQuadratureOrder =
      parameterList.template get<int>("options.quadrature.near.doubleOrder");
  const int testQuadOrder = nearDoubleQuadratureOrder;
  const int trialQuadOrder = nearDoubleQuadratureOrder;
  std::cout << "testQuadOrder = " << testQuadOrder << ", " << std::flush;
  std::cout << "trialQuadOrder = " << trialQuadOrder << std::endl;

  Matrix<CoordinateType> localTestQuadPoints, localTrialQuadPoints;
  std::vector<CoordinateType> testQuadWeights, trialQuadWeights;

  Fiber::fillSingleQuadraturePointsAndWeights(
      testGrid->idxCount(), testQuadOrder, localTestQuadPoints, testQuadWeights);
  Fiber::fillSingleQuadraturePointsAndWeights(
      trialGrid->idxCount(), trialQuadOrder, localTrialQuadPoints, trialQuadWeights);

  const unsigned int testPointCount = localTestQuadPoints.cols();
  const unsigned int trialPointCount = localTrialQuadPoints.cols();

  const unsigned int maxActiveElemPairCount =
      getMaxActiveElemPairCount<ResultType>(
          testGrid, trialGrid,
          testDofCount, trialDofCount,
          testPointCount, trialPointCount,
          cudaOptions.isElemDataCachingEnabled());
  std::cout << "maxActiveElemPairCount = " << maxActiveElemPairCount << ", " << std::flush;

  const unsigned int chunkCount = std::max(
      static_cast<int>((elemPairCount-1)/maxActiveElemPairCount+1), 1);
  std::cout << "chunkCount = " << chunkCount << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  thrust::host_vector<ResultType> h_regularResultEven(
      maxActiveElemPairCount * testDofCount * trialDofCount);
  thrust::host_vector<ResultType> h_regularResultOdd(
      maxActiveElemPairCount * testDofCount * trialDofCount);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for host vector allocation = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  Fiber::CudaIntegrator<
      BasisFunctionType, KernelType, ResultType> cudaIntegrator(
          localTestQuadPoints, localTrialQuadPoints,
          testQuadWeights, trialQuadWeights,
          testShapeset, trialShapeset,
          testGrid, trialGrid, *kernels,
          cudaOptions.isElemDataCachingEnabled(), cudaOptions.streamCount());

  bool singularIntegralCachingEnabled =
      context.assemblyOptions().isSingularIntegralCachingEnabled();
  if (!singularIntegralCachingEnabled)
    throw std::invalid_argument(
        "CudaDenseGlobalAssembler::assembleDetachedWeakForm(): "
        "singular integral caching has to be enabled");

  // Get reference to singular integral cache
  const Cache &cache = assembler.cache();

  // Loop over element pair chunks
  for (int chunk = 0; chunk < chunkCount; ++chunk) {

    std::cout << "chunk = " << chunk << std::endl;

    thrust::host_vector<ResultType>* h_regularResult;
    if (chunk % 2 == 0)
      h_regularResult = &h_regularResultEven;
    else
      h_regularResult = &h_regularResultOdd;

    unsigned int chunkElemPairCount;
    if (chunk == chunkCount-1) {
      chunkElemPairCount = testElemPairIndices.end()
          - (testElemPairIndices.begin() + chunk * maxActiveElemPairCount);
    } else {
      chunkElemPairCount = maxActiveElemPairCount;
    }
    std::cout << "chunkElemPairCount = " << chunkElemPairCount << std::endl;

    taskGroupGlobal.run([chunk, chunkCount, maxActiveElemPairCount,
         chunkElemPairCount, testDofCount, trialDofCount,
         &testElemPairIndices, &trialElemPairIndices,
         &cudaIntegrator, h_regularResult]{

      std::vector<int>::iterator startElemPairTestIndex =
          testElemPairIndices.begin() + chunk * maxActiveElemPairCount;
      std::vector<int>::iterator startElemPairTrialIndex =
          trialElemPairIndices.begin() + chunk * maxActiveElemPairCount;

      std::vector<int>::iterator endElemPairTestIndex;
      std::vector<int>::iterator endElemPairTrialIndex;
      if (chunk == chunkCount-1) {
        endElemPairTestIndex = testElemPairIndices.end();
        endElemPairTrialIndex = trialElemPairIndices.end();
      } else {
        endElemPairTestIndex = testElemPairIndices.begin()
            + (chunk+1) * maxActiveElemPairCount;
        endElemPairTrialIndex = trialElemPairIndices.begin()
            + (chunk+1) * maxActiveElemPairCount;
      }

      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

      // Evaluate regular integrals over selected element pairs
      cudaIntegrator.integrate(startElemPairTestIndex, endElemPairTestIndex,
                               startElemPairTrialIndex, endElemPairTrialIndex,
                               h_regularResult->begin());

      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      std::cout << "Time for CudaIntegrator::integrate() = "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                << " ms" << std::endl;
    });
    taskGroupGlobal.wait();

    taskGroupGlobal.run([chunk, chunkCount, maxActiveElemPairCount,
        chunkElemPairCount, testDofCount, trialDofCount,
        &testElemPairIndices, &trialElemPairIndices,
        &testGlobalDofs, &trialGlobalDofs,
        &testLocalDofWeights, &trialLocalDofWeights,
        h_regularResult, cache, &result, &mutex]{

      const int offset = chunk * maxActiveElemPairCount;

      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

      // Global assembly
      tbb::spin_mutex::scoped_lock lock(mutex);

      for (int chunkElemPair = 0; chunkElemPair < chunkElemPairCount; ++chunkElemPair) {

        const int testIndex = testElemPairIndices[offset + chunkElemPair];
        const int trialIndex = trialElemPairIndices[offset + chunkElemPair];

        // Try to find matrix in singular integral cache
        const Matrix<ResultType> *cachedLocalWeakForm = 0;
        for (size_t n = 0; n < cache.extent(0); ++n)
          if (cache(n, trialIndex).first == testIndex) {
            cachedLocalWeakForm = &cache(n, trialIndex).second;
            break;
          }

        const int testDofCount = testGlobalDofs[testIndex].size();
        const int trialDofCount = trialGlobalDofs[trialIndex].size();

        // Add the integrals to appropriate entries in the operator's matrix
        for (int testDof = 0; testDof < testDofCount; ++testDof) {

          int testGlobalDof = testGlobalDofs[testIndex][testDof];
          if (testGlobalDof < 0)
            continue;

          for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {

            int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
            if (trialGlobalDof < 0)
              continue;

            assert(std::abs(testLocalDofWeights[testIndex][testDof]) > 0.);
            assert(std::abs(trialLocalDofWeights[trialIndex][trialDof]) > 0.);

            if (cachedLocalWeakForm) { // Matrix found in cache
              result(testGlobalDof, trialGlobalDof) +=
                conj(testLocalDofWeights[testIndex][testDof]) *
                trialLocalDofWeights[trialIndex][trialDof] *
                (*cachedLocalWeakForm)(testDof, trialDof);
            } else {
              result(testGlobalDof, trialGlobalDof) +=
                conj(testLocalDofWeights[testIndex][testDof]) *
                trialLocalDofWeights[trialIndex][trialDof] *
                (*h_regularResult)[chunkElemPair * testDofCount * trialDofCount
                                + testDof * trialDofCount
                                + trialDof];
            }
          }
        }
      }
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      std::cout << "Time for regular global result assembly = "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                << " ms" << std::endl;
    });
  }
  taskGroupGlobal.wait();

  cudaProfilerStop();

  if (result.rows() < 10 && result.cols() < 10) {
    std::cout << "result (cudadense) = " << std::endl;
    std::cout << result << std::endl;
  }

  // Create and return a discrete operator represented by the matrix that
  // has just been calculated
  return std::unique_ptr<DiscreteBoundaryOperator<ResultType>>(
      new DiscreteDenseBoundaryOperator<ResultType>(result));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaDenseGlobalAssembler);

} // namespace Bempp
