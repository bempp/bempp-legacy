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

#include "../common/types.hpp"
#include "../common/complex_aux.hpp"

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

template <typename BasisFunctionType>
class GetSortedElementPairsLoopBody {
public:

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
//  typedef tbb::spin_mutex MutexType;
  typedef tbb::mutex MutexType;

  typedef std::pair<Matrix<CoordinateType>, std::vector<CoordinateType>> QuadDataset;
  typedef std::tuple<const QuadDataset*, const QuadDataset*, const bool> Integrator;
  typedef tbb::concurrent_unordered_map<
      Fiber::DoubleQuadratureDescriptor, Integrator*> IntegratorMap;

  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;
  typedef std::pair<const Shapeset*, const Shapeset*> ShapesetPair;

  typedef std::pair<const Integrator*, ShapesetPair> QuadVariant;

  GetSortedElementPairsLoopBody(
      const std::vector<int> &testIndices, const std::vector<int> &trialIndices,
      std::vector<const Shapeset*> &testShapesets,
      std::vector<const Shapeset*> &trialShapesets,
      shared_ptr<const Fiber::QuadratureDescriptorSelectorForIntegralOperators<
      CoordinateType>> quadDescSelector,
      shared_ptr<const Fiber::DoubleQuadratureRuleFamily<CoordinateType>>
      quadRuleFamily,
      tbb::concurrent_unordered_map<
          Fiber::DoubleQuadratureDescriptor, Integrator*> integrators,
      const CoordinateType nominalDistance,
      std::vector<std::vector<QuadVariant>> &quadVariants,
      MutexType &mutex)
      : m_testIndices(testIndices), m_trialIndices(trialIndices),
        m_testShapesets(testShapesets), m_trialShapesets(trialShapesets),
        m_quadDescSelector(quadDescSelector), m_quadRuleFamily(quadRuleFamily),
        m_integrators(integrators), m_nominalDistance(nominalDistance),
        m_quadVariants(quadVariants), m_mutex(mutex) {}

  void operator()(const tbb::blocked_range<int> &r) const {

    const int trialIndexCount = m_trialIndices.size();
    for (int i = r.begin(); i != r.end(); ++i) {
      for (int j = 0; j < trialIndexCount; ++j) {

        const int testIndex = m_testIndices[i];
        const int trialIndex = m_trialIndices[j];

        const Integrator *integrator = &selectIntegrator(testIndex, trialIndex);

        const Shapeset* testShapeset = m_testShapesets[testIndex];
        const Shapeset* trialShapeset = m_trialShapesets[trialIndex];

        ShapesetPair* shapesetPair = new ShapesetPair(testShapeset, trialShapeset);
        QuadVariant* quadVariant = new QuadVariant(integrator, *shapesetPair);
        m_quadVariants[i][j] = *quadVariant;
      }
    }
  }

private:

  const Integrator& selectIntegrator(
      const int testIndex, const int trialIndex) const {

    const Fiber::DoubleQuadratureDescriptor desc =
        m_quadDescSelector->quadratureDescriptor(
            testIndex, trialIndex, m_nominalDistance);
    return getIntegrator(desc);
  }

  // TODO: creates memory leak
  const Integrator& getIntegrator(
      const Fiber::DoubleQuadratureDescriptor &desc) const {

    typename IntegratorMap::iterator it = m_integrators.find(desc);
    if (it == m_integrators.end()) {

//      tbb::mutex::scoped_lock lock(m_integratorCreationMutex);
      MutexType::scoped_lock lock(m_mutex);
      it = m_integrators.find(desc);

      if (it == m_integrators.end()) {

        // Integrator doesn't exist yet and must be created
        Matrix<CoordinateType> testPoints, trialPoints;
        std::vector<CoordinateType> testWeights, trialWeights;
        bool isTensor;
        m_quadRuleFamily->fillQuadraturePointsAndWeights(
            desc, testPoints, trialPoints, testWeights, trialWeights, isTensor);

        const QuadDataset *testQuadDataset = new QuadDataset(testPoints, testWeights);
        const QuadDataset *trialQuadDataset = new QuadDataset(trialPoints, trialWeights);
        Integrator *integrator = new Integrator(
            testQuadDataset, trialQuadDataset,
            isTensor);

        // Attempt to insert the newly created integrator into the map
        std::pair<typename IntegratorMap::iterator, bool> result =
            m_integrators.insert(std::make_pair(desc, integrator));

        if (result.second)
          // Insertion succeeded. The newly created integrator will be deleted in
          // our own destructor
          ;
        else
          // Insertion failed -- another thread was faster. Delete the newly
          // created integrator
          delete integrator;

        // Return pointer to the integrator that ended up in the map
        it = result.first;
      }
    }
    return *it->second;
  }

  const std::vector<int> &m_testIndices;
  const std::vector<int> &m_trialIndices;

  std::vector<const Shapeset*> &m_testShapesets;
  std::vector<const Shapeset*> &m_trialShapesets;

  shared_ptr<const Fiber::QuadratureDescriptorSelectorForIntegralOperators<
  CoordinateType>> &m_quadDescSelector;
  shared_ptr<const Fiber::DoubleQuadratureRuleFamily<CoordinateType>>
  &m_quadRuleFamily;

  IntegratorMap &m_integrators;

  const CoordinateType m_nominalDistance;

  std::vector<std::vector<QuadVariant>> &m_quadVariants;

  // mutex must be mutable because we need to lock and unlock it
  MutexType &m_mutex;
};

// Get sorted element pair index vectors regarding regular and singular integrals
// as well as quadrature order and shapeset combinations for regular ones
template <typename BasisFunctionType>
void getSortedElementPairs(
    const Space<BasisFunctionType> &testSpace,
    const Space<BasisFunctionType> &trialSpace,
    std::vector<int> &testIndices,
    std::vector<int> &trialIndices,
    shared_ptr<const Fiber::QuadratureDescriptorSelectorForIntegralOperators<
    typename ScalarTraits<BasisFunctionType>::RealType>> quadDescSelector,
    shared_ptr<const Fiber::DoubleQuadratureRuleFamily<
    typename ScalarTraits<BasisFunctionType>::RealType>> quadRuleFamily,
    std::vector<std::vector<int>> &regularElemPairTestIndices,
    std::vector<std::vector<int>> &regularElemPairTrialIndices,
    std::vector<std::pair<const std::tuple<const std::pair<Matrix<typename ScalarTraits<BasisFunctionType>::RealType>, std::vector<typename ScalarTraits<BasisFunctionType>::RealType>>*, const std::pair<Matrix<typename ScalarTraits<BasisFunctionType>::RealType>, std::vector<typename ScalarTraits<BasisFunctionType>::RealType>>*, const bool>*, std::pair<const Fiber::Shapeset<BasisFunctionType>*, const Fiber::Shapeset<BasisFunctionType>*>>>
    &regularQuadVariants,
    std::vector<int> &singularElemPairTestIndices,
    std::vector<int> &singularElemPairTrialIndices,
    typename ScalarTraits<BasisFunctionType>::RealType nominalDistance = -1.) {

  const int testIndexCount = testIndices.size();
  const int trialIndexCount = trialIndices.size();

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;
  typedef std::pair<const Shapeset*, const Shapeset*> ShapesetPair;
  std::vector<const Shapeset*> testShapesets, trialShapesets;
  getAllShapesets(testSpace, testShapesets);
  if (&testSpace == &trialSpace) {
    trialShapesets = testShapesets;
  } else {
    getAllShapesets(trialSpace, trialShapesets);
  }

  typedef std::pair<Matrix<CoordinateType>, std::vector<CoordinateType>> QuadDataset;
  typedef std::tuple<const QuadDataset*, const QuadDataset*, const bool> Integrator;

  typedef std::pair<const Integrator*, ShapesetPair> QuadVariant;
  std::vector<std::vector<QuadVariant>> quadVariants(testIndexCount);
  for (int i = 0; i < testIndexCount; ++i)
    quadVariants[i].resize(trialIndexCount);

  typedef tbb::concurrent_unordered_map<
      Fiber::DoubleQuadratureDescriptor, Integrator*> IntegratorMap;
  IntegratorMap integrators;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  typedef GetSortedElementPairsLoopBody<BasisFunctionType> Body;
  typename Body::MutexType mutex;
  {
    Fiber::SerialBlasRegion region;
    tbb::parallel_for(tbb::blocked_range<int>(0, testIndexCount),
                      Body(testIndices, trialIndices,
                           testShapesets, trialShapesets,
                           quadDescSelector, quadRuleFamily, integrators,
                           nominalDistance, quadVariants, mutex));
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for GetSortedElementPairsLoopBody() = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  // Integration will proceed in batches of element pairs having the same
  // "quadrature variant", i.e. quadrature datasets and shapesets

  begin = std::chrono::steady_clock::now();

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants;
  for (int i = 0; i < quadVariants.size(); ++i)
    uniqueQuadVariants.insert(quadVariants[i].begin(), quadVariants[i].end());

  end = std::chrono::steady_clock::now();
  std::cout << "Time for creating set of unique quad variants = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  begin = std::chrono::steady_clock::now();

//  std::vector<tbb::concurrent_vector<int>> tbb_regularElemPairTestIndices;
//  std::vector<tbb::concurrent_vector<int>> tbb_regularElemPairTrialIndices;
//  tbb::concurrent_vector<int> tbb_singularElemPairTestIndices;
//  tbb::concurrent_vector<int> tbb_singularElemPairTrialIndices;
//  tbb::queuing_mutex pushBackMutex;

  // Now loop over unique quadrature variants
  unsigned int regularQuadVariantCount = 0;
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {

    const QuadVariant activeQuadVariant = *it;
    const Integrator &activeIntegrator = *it->first;
    ShapesetPair activeBasisPair = it->second;
    const bool isRegular = std::get<2>(activeIntegrator);

    // Find all the element pairs which should be considered according to the
    // current quadrature variant
    if (isRegular) {

      regularQuadVariants.push_back(activeQuadVariant);
      regularElemPairTestIndices.resize(regularQuadVariantCount+1);
      regularElemPairTrialIndices.resize(regularQuadVariantCount+1);

      for (int i = 0; i < testIndexCount; ++i)
//      tbb::parallel_for(tbb::blocked_range<int>(0, testIndexCount),
//          [trialIndexCount, regularQuadVariantCount,
//           &testIndices, &trialIndices,
//           &quadVariants, &activeQuadVariant,
//           &tbb_regularElemPairTestIndices, &tbb_regularElemPairTrialIndices,
//           &pushBackMutex](const tbb::blocked_range<int>& r) {
//
//        for (int i = r.begin(); i != r.end(); ++i) {
          for (int j = 0; j < trialIndexCount; ++j) {
            if (quadVariants[i][j] == activeQuadVariant) {

              // Element pair related to regular integrals
              regularElemPairTestIndices[regularQuadVariantCount].push_back(testIndices[i]);
              regularElemPairTrialIndices[regularQuadVariantCount].push_back(trialIndices[j]);
            }
          }
//        }
//      });
      regularQuadVariantCount++;

    } else {

      for (int i = 0; i < testIndexCount; ++i)
//      tbb::parallel_for(tbb::blocked_range<int>(0, testIndexCount),
//          [trialIndexCount,
//           &testIndices, &trialIndices,
//           &quadVariants, &activeQuadVariant,
//           &tbb_singularElemPairTestIndices, &tbb_singularElemPairTrialIndices,
//           &pushBackMutex](const tbb::blocked_range<int>& r) {
//
//        for (int i = r.begin(); i != r.end(); ++i) {
          for (int j = 0; j < trialIndexCount; ++j) {
            if (quadVariants[i][j] == activeQuadVariant) {

              // Element pair related to singular integrals
              singularElemPairTestIndices.push_back(testIndices[i]);
              singularElemPairTrialIndices.push_back(trialIndices[j]);
            }
          }
//        }
//      });
    }
  }
//  regularElemPairTestIndices.resize(regularQuadVariantCount);
//  regularElemPairTrialIndices.resize(regularQuadVariantCount);
//  for (int regularQuadVariant = 0; regularQuadVariant < regularQuadVariantCount; ++regularQuadVariant) {
//    regularElemPairTestIndices[regularQuadVariant].resize(
//        tbb_regularElemPairTestIndices[regularQuadVariant].size());
//    regularElemPairTrialIndices[regularQuadVariant].resize(
//        tbb_regularElemPairTrialIndices[regularQuadVariant].size());
//    regularElemPairTestIndices[regularQuadVariant].assign(
//        tbb_regularElemPairTestIndices[regularQuadVariant].begin(),
//        tbb_regularElemPairTestIndices[regularQuadVariant].end());
//    regularElemPairTrialIndices[regularQuadVariant].assign(
//        tbb_regularElemPairTrialIndices[regularQuadVariant].begin(),
//        tbb_regularElemPairTrialIndices[regularQuadVariant].end());
//  }
//  singularElemPairTestIndices.resize(tbb_singularElemPairTestIndices.size());
//  singularElemPairTrialIndices.resize(tbb_singularElemPairTrialIndices.size());
//  singularElemPairTestIndices.assign(tbb_singularElemPairTestIndices.begin(),
//      tbb_singularElemPairTestIndices.end());
//  singularElemPairTrialIndices.assign(tbb_singularElemPairTrialIndices.begin(),
//      tbb_singularElemPairTrialIndices.end());

  end = std::chrono::steady_clock::now();
  std::cout << "Time for looping over unique quad variants = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;
}

// Determine the maximum number of element pairs in a chunk treated on the
// device according to hardware specs
template <typename ResultType>
unsigned int getMaxActiveElemPairCount(shared_ptr<CudaGrid> testGrid,
                                       shared_ptr<CudaGrid> trialGrid,
                                       const unsigned int testDofCount,
                                       const unsigned int trialDofCount,
                                       const unsigned int testPointCount,
                                       const unsigned int trialPointCount,
                                       size_t reserveMem = 1e09,
                                       bool cacheElemData = false) {

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

  // Get element pair index vectors and regular integration data
  std::vector<std::vector<int>> regularElemPairTestIndices;
  std::vector<std::vector<int>> regularElemPairTrialIndices;
  std::vector<QuadVariant> regularQuadVariants;
  std::vector<int> singularElemPairTestIndices;
  std::vector<int> singularElemPairTrialIndices;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  getSortedElementPairs(
      testSpace, trialSpace, testIndices, trialIndices,
      assembler.quadDescSelector(), assembler.quadRuleFamily(),
      regularElemPairTestIndices, regularElemPairTrialIndices,
      regularQuadVariants,
      singularElemPairTestIndices, singularElemPairTrialIndices);

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for getSortedElementPairs() = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;

  std::cout << "number of singular element pairs: " << singularElemPairTestIndices.size() << std::endl;
  std::cout << "number of regular element pairs: " << std::flush;
  for (int i = 0; i < regularElemPairTestIndices.size(); ++i)
    std::cout << regularElemPairTestIndices[i].size() << " " << std::flush;
  std::cout << std::endl;

//  std::cout << "singularElemPairIndices = " << std::endl;
//  for (int i = 0; i < singularElemPairTestIndices.size(); ++i) {
//    std::cout << singularElemPairTestIndices[i] << " "
//              << singularElemPairTrialIndices[i] << std::endl;
//  }
//  std::cout << std::endl;
//
//  std::cout << "regularElemPairIndices = " << std::endl;
//  for (int i = 0; i < regularElemPairTestIndices[0].size(); ++i) {
//    std::cout << regularElemPairTestIndices[0][i] << " "
//              << regularElemPairTrialIndices[0][i] << std::endl;
//  }
//  std::cout << std::endl;

  std::cout << "number of quadrature points = " << std::endl;
  for (int i = 0; i < regularQuadVariants.size(); ++i) {
    std::cout << (std::get<0>(*(regularQuadVariants[i].first))->second).size() << " " << std::flush;
  }
  std::cout << std::endl;
  for (int i = 0; i < regularQuadVariants.size(); ++i) {
    std::cout << (std::get<1>(*(regularQuadVariants[i].first))->second).size() << " " << std::flush;
  }
  std::cout << std::endl;

  std::cout << "number of shape functions: " << std::endl;
  for (int i = 0; i < regularQuadVariants.size(); ++i) {
    std::cout << regularQuadVariants[i].second.first->size() << " " << std::flush;
  }
  std::cout << std::endl;
  for (int i = 0; i < regularQuadVariants.size(); ++i) {
    std::cout << regularQuadVariants[i].second.second->size() << " " << std::flush;
  }
  std::cout << std::endl;

  tbb::task_group taskGroupGlobal;
  // TODO: Most suitable mutex type?
  tbb::spin_mutex mutex;

  taskGroupGlobal.run([&singularElemPairTestIndices, &singularElemPairTrialIndices,
      &testGlobalDofs, &trialGlobalDofs,
      &testLocalDofWeights, &trialLocalDofWeights,
      &assembler, &result, &mutex]{

    std::vector<Matrix<ResultType>> singularResult;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Evaluate singular integrals over selected element pairs
    assembler.evaluateLocalWeakForms(singularElemPairTestIndices,
                                     singularElemPairTrialIndices,
                                     singularResult);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time for evaluateSingularWeakForms() = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << " ms" << std::endl;

    begin = std::chrono::steady_clock::now();

    // Global assembly
    tbb::spin_mutex::scoped_lock lock(mutex);

    // Loop over singular element pairs
    for (int singularElemPair = 0; singularElemPair < singularElemPairTestIndices.size(); ++singularElemPair) {

      const int testIndex = singularElemPairTestIndices[singularElemPair];
      const int trialIndex = singularElemPairTrialIndices[singularElemPair];

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

          result(testGlobalDof, trialGlobalDof) +=
            conj(testLocalDofWeights[testIndex][testDof]) *
            trialLocalDofWeights[trialIndex][trialDof] *
            singularResult[singularElemPair](testDof, trialDof);
        }
      }
    }

    end = std::chrono::steady_clock::now();
    std::cout << "Time for singular global assembly = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << " ms" << std::endl;
  });

  taskGroupGlobal.run([&regularElemPairTestIndices, &regularElemPairTrialIndices,
      &testGlobalDofs, &trialGlobalDofs,
      &testLocalDofWeights, &trialLocalDofWeights,
      &testSpace, &trialSpace,
      &regularQuadVariants, &assembler, &result, &mutex]{

    // TODO: How to determine the correct kernel type?
    typedef CoordinateType KernelType;         // Real kernel
//        typedef ResultType KernelType;             // Complex kernel
    shared_ptr<const Fiber::CollectionOfKernels<KernelType>> kernels;
//        assembler.getKernels(kernels);

    cudaProfilerStart();

    // Push raw grid data to the device
    shared_ptr<CudaGrid> testGrid = testSpace.grid()->pushToDevice(0);
    shared_ptr<CudaGrid> trialGrid = trialSpace.grid()->pushToDevice(0);

    const unsigned int quadVariantCount = regularQuadVariants.size();

    std::vector<unsigned int> elemPairCounts(quadVariantCount, 0);
    std::vector<unsigned int> maxActiveElemPairCounts(quadVariantCount, 0);
    std::vector<unsigned int> chunkCounts(quadVariantCount, 0);

    unsigned int maxTestDofCount = 0, maxTrialDofCount = 0;

    for (int quadVariant = 0; quadVariant < quadVariantCount; ++quadVariant) {

      elemPairCounts[quadVariant] =
          regularElemPairTestIndices[quadVariant].size();

      const unsigned int testPointCount =
          (std::get<0>(*(regularQuadVariants[quadVariant].first))->second).size();
      const unsigned int trialPointCount =
          (std::get<1>(*(regularQuadVariants[quadVariant].first))->second).size();

      const unsigned int testDofCount =
          (*(regularQuadVariants[quadVariant].second.first)).size();
      const unsigned int trialDofCount =
          (*(regularQuadVariants[quadVariant].second.second)).size();

      if (testDofCount > maxTestDofCount) maxTestDofCount = testDofCount;
      if (trialDofCount > maxTrialDofCount) maxTrialDofCount = trialDofCount;

      maxActiveElemPairCounts[quadVariant] =
          getMaxActiveElemPairCount<ResultType>(
              testGrid, trialGrid,
              testDofCount, trialDofCount,
              testPointCount, trialPointCount);
//      maxActiveElemPairCounts[0] = 3;

      chunkCounts[quadVariant] = std::max(
          static_cast<int>((elemPairCounts[quadVariant]-1)
              / maxActiveElemPairCounts[quadVariant]+1), 1);
    }
    thrust::host_vector<ResultType> h_regularResultEven(
        *std::max_element(maxActiveElemPairCounts.begin(),
            maxActiveElemPairCounts.end()) * maxTestDofCount * maxTrialDofCount);
    thrust::host_vector<ResultType> h_regularResultOdd(
        *std::max_element(maxActiveElemPairCounts.begin(),
            maxActiveElemPairCounts.end()) * maxTestDofCount * maxTrialDofCount);

    unsigned int maxActiveElemPairCount = 0;
    unsigned int localChunkCount = 0;
    unsigned int testDofCount = 0;
    unsigned int trialDofCount = 0;

    shared_ptr<Fiber::CudaIntegrator<
        BasisFunctionType, KernelType, ResultType>> cudaIntegrator;

    const unsigned int totalChunkCount =
        std::accumulate(chunkCounts.begin(), chunkCounts.end(), 0);
    std::cout << "totalChunkCount = " << totalChunkCount << std::endl;

    tbb::task_group taskGroupRegular;

    int quadVariant = 0, localChunk = 0;
    for (int chunk = 0; chunk < totalChunkCount; ++chunk) {

      std::cout << "chunk = " << chunk << std::endl;
      std::cout << "quadVariant = " << quadVariant << std::endl;
      std::cout << "localChunk = " << localChunk << std::endl;

      thrust::host_vector<ResultType>* h_regularResult;
      if (chunk % 2 == 0)
        h_regularResult = &h_regularResultEven;
      else
        h_regularResult = &h_regularResultOdd;

      // First time this quadrature variant is used; create new integrator
      if (localChunk == 0) {

        Matrix<CoordinateType> localTestQuadPoints =
            std::get<0>(*(regularQuadVariants[quadVariant].first))->first;
        Matrix<CoordinateType> localTrialQuadPoints =
            std::get<1>(*(regularQuadVariants[quadVariant].first))->first;
        std::vector<CoordinateType> testQuadWeights =
            std::get<0>(*(regularQuadVariants[quadVariant].first))->second;
        std::vector<CoordinateType> trialQuadWeights =
            std::get<1>(*(regularQuadVariants[quadVariant].first))->second;

        const Shapeset &testShapeset =
            *(regularQuadVariants[quadVariant].second.first);
        const Shapeset &trialShapeset =
            *(regularQuadVariants[quadVariant].second.second);

        testDofCount = testShapeset.size();
        trialDofCount = trialShapeset.size();

        // Create chunks of element pairs according to a maximum number of
        // element pairs active on the device
        maxActiveElemPairCount = maxActiveElemPairCounts[quadVariant];
        std::cout << "maxActiveElemPairCount = " << maxActiveElemPairCount << std::endl;
        std::cout << "localElemPairCount = " << elemPairCounts[quadVariant] << std::endl;

        cudaIntegrator = boost::make_shared<Fiber::CudaIntegrator<
            BasisFunctionType, KernelType, ResultType>>(
                localTestQuadPoints, localTrialQuadPoints,
                testQuadWeights, trialQuadWeights,
                testShapeset, trialShapeset,
                testGrid, trialGrid,
                *kernels);

        localChunkCount = chunkCounts[quadVariant];
        std::cout << "localChunkCount = " << localChunkCount << std::endl;
      }

      taskGroupRegular.run([quadVariant, localChunk, localChunkCount,
           maxActiveElemPairCount, testDofCount, trialDofCount,
           &regularElemPairTestIndices, &regularElemPairTrialIndices,
           &cudaIntegrator, h_regularResult]{

        std::vector<int>::iterator startElemPairTestIndex =
            regularElemPairTestIndices[quadVariant].begin()
            + localChunk * maxActiveElemPairCount;
        std::vector<int>::iterator startElemPairTrialIndex =
            regularElemPairTrialIndices[quadVariant].begin()
            + localChunk * maxActiveElemPairCount;

        std::vector<int>::iterator endElemPairTestIndex;
        std::vector<int>::iterator endElemPairTrialIndex;
        if (localChunk == localChunkCount-1) {
          endElemPairTestIndex =
              regularElemPairTestIndices[quadVariant].end();
          endElemPairTrialIndex =
              regularElemPairTrialIndices[quadVariant].end();
        } else {
          endElemPairTestIndex =
              regularElemPairTestIndices[quadVariant].begin()
              + (localChunk+1) * maxActiveElemPairCount;
          endElemPairTrialIndex =
              regularElemPairTrialIndices[quadVariant].begin()
              + (localChunk+1) * maxActiveElemPairCount;
        }

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        // Evaluate regular integrals over selected element pairs
        cudaIntegrator->integrate(startElemPairTestIndex, endElemPairTestIndex,
                                  startElemPairTrialIndex, endElemPairTrialIndex,
                                  h_regularResult->begin());

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time for CudaIntegrator::integrate() = "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                  << " ms" << std::endl;
      });
      taskGroupRegular.wait();

      taskGroupRegular.run([quadVariant, localChunk, localChunkCount,
          maxActiveElemPairCount, testDofCount, trialDofCount,
          &regularElemPairTestIndices, &regularElemPairTrialIndices,
          &testGlobalDofs, &trialGlobalDofs,
          &testLocalDofWeights, &trialLocalDofWeights,
          h_regularResult, &result, &mutex]{

          const int offset = localChunk * maxActiveElemPairCount;
          unsigned int chunkElemPairCount;
          if (localChunk == localChunkCount-1) {
            chunkElemPairCount = regularElemPairTestIndices[quadVariant].end()
                - (regularElemPairTestIndices[quadVariant].begin()
                   + localChunk * maxActiveElemPairCount);
          } else {
            chunkElemPairCount = maxActiveElemPairCount;
          }

          std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

          // Global assembly
          tbb::spin_mutex::scoped_lock lock(mutex);

          for (int chunkElemPair = 0; chunkElemPair < chunkElemPairCount; ++chunkElemPair) {

            const int testIndex =
                regularElemPairTestIndices[quadVariant][offset + chunkElemPair];
            const int trialIndex =
                regularElemPairTrialIndices[quadVariant][offset + chunkElemPair];

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

                result(testGlobalDof, trialGlobalDof) +=
                  conj(testLocalDofWeights[testIndex][testDof]) *
                  trialLocalDofWeights[trialIndex][trialDof] *
                  (*h_regularResult)[chunkElemPair * testDofCount * trialDofCount
                                  + testDof * trialDofCount
                                  + trialDof];
              }
            }
          }
          std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
          std::cout << "Time for regular global result assembly = "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
                    << " ms" << std::endl;
        });
      if (localChunk == localChunkCount-1) {
        localChunk = 0;
        quadVariant++;
      } else {
        localChunk++;
      }
    }
    taskGroupRegular.wait();

    cudaProfilerStop();
  });

  taskGroupGlobal.wait();

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
