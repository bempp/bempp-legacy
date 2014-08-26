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

#include "../common/common.hpp"

// Keep IDEs happy
#include "default_local_assembler_for_integral_operators_on_surfaces.hpp"

#include "double_quadrature_rule_family.hpp"
#include "nonseparable_numerical_test_kernel_trial_integrator.hpp"
#include "quadrature_descriptor_selector_for_integral_operators.hpp"
#include "separable_numerical_test_kernel_trial_integrator.hpp"
#include "serial_blas_region.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include "../common/auto_timer.hpp"

namespace Fiber {

namespace {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class SingularIntegralCalculatorLoopBody {
public:
  typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
  Integrator;
  typedef typename Integrator::ElementIndexPair ElementIndexPair;

  SingularIntegralCalculatorLoopBody(
      const Integrator &activeIntegrator,
      const std::vector<ElementIndexPair> &activeElementPairs,
      const Shapeset<BasisFunctionType> &activeTestShapeset,
      const Shapeset<BasisFunctionType> &activeTrialShapeset,
      const std::vector<arma::Mat<ResultType> *> &localResult)
      : m_activeIntegrator(activeIntegrator),
        m_activeElementPairs(activeElementPairs),
        m_activeTestBasis(activeTestShapeset),
        m_activeTrialBasis(activeTrialShapeset), m_localResult(localResult) {}

  void operator()(const tbb::blocked_range<size_t> &r) const {
    // copy the relevant subset of m_activeElementPairs into
    // localActiveElementPairs
    std::vector<ElementIndexPair> localActiveElementPairs(
        &m_activeElementPairs[r.begin()], &m_activeElementPairs[r.end()]);
    std::vector<arma::Mat<ResultType> *> localLocalResult(
        &m_localResult[r.begin()], &m_localResult[r.end()]);
    m_activeIntegrator.integrate(localActiveElementPairs, m_activeTestBasis,
                                 m_activeTrialBasis, localLocalResult);
  }

private:
  const Integrator &m_activeIntegrator;
  const std::vector<ElementIndexPair> &m_activeElementPairs;
  const Shapeset<BasisFunctionType> &m_activeTestBasis;
  const Shapeset<BasisFunctionType> &m_activeTrialBasis;
  const std::vector<arma::Mat<ResultType> *> &m_localResult;
};

} // namespace

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    DefaultLocalAssemblerForIntegralOperatorsOnSurfaces(
        const shared_ptr<const GeometryFactory> &testGeometryFactory,
        const shared_ptr<const GeometryFactory> &trialGeometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType>> &
            testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType>> &
            trialRawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType> *>> &testShapesets,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType> *>> &trialShapesets,
        const shared_ptr<const CollectionOfShapesetTransformations<
            CoordinateType>> &testTransformations,
        const shared_ptr<const CollectionOfKernels<KernelType>> &kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<
            CoordinateType>> &trialTransformations,
        const shared_ptr<const TestKernelTrialIntegral<
            BasisFunctionType, KernelType, ResultType>> &integral,
        const shared_ptr<const OpenClHandler> &openClHandler,
        const ParallelizationOptions &parallelizationOptions,
        VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals,
        const shared_ptr<const QuadratureDescriptorSelectorForIntegralOperators<
            CoordinateType>> &quadDescSelector,
        const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
            quadRuleFamily)
    : m_testGeometryFactory(testGeometryFactory),
      m_trialGeometryFactory(trialGeometryFactory),
      m_testRawGeometry(testRawGeometry), m_trialRawGeometry(trialRawGeometry),
      m_testShapesets(testShapesets), m_trialShapesets(trialShapesets),
      m_testTransformations(testTransformations), m_kernels(kernels),
      m_trialTransformations(trialTransformations), m_integral(integral),
      m_openClHandler(openClHandler),
      m_parallelizationOptions(parallelizationOptions),
      m_verbosityLevel(verbosityLevel), m_quadDescSelector(quadDescSelector),
      m_quadRuleFamily(quadRuleFamily) {
  Utilities::checkConsistencyOfGeometryAndShapesets(*testRawGeometry,
                                                    *testShapesets);
  Utilities::checkConsistencyOfGeometryAndShapesets(*trialRawGeometry,
                                                    *trialShapesets);

  if (cacheSingularIntegrals)
    cacheSingularLocalWeakForms();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::~DefaultLocalAssemblerForIntegralOperatorsOnSurfaces() {
  // Note: obviously the destructor is assumed to be called only after
  // all threads have ceased using the assembler!

  for (typename IntegratorMap::const_iterator it =
           m_testKernelTrialIntegrators.begin();
       it != m_testKernelTrialIntegrators.end(); ++it)
    delete it->second;
  m_testKernelTrialIntegrators.clear();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
inline bool DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::testAndTrialGridsAreIdentical() const {
  return m_testRawGeometry.get() == m_trialRawGeometry.get();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    evaluateLocalWeakForms(CallVariant callVariant,
                           const std::vector<int> &elementIndicesA,
                           int elementIndexB, LocalDofIndex localDofIndexB,
                           std::vector<arma::Mat<ResultType>> &result,
                           CoordinateType nominalDistance) {
  typedef Shapeset<BasisFunctionType> Shapeset;

  const int elementACount = elementIndicesA.size();
  result.resize(elementACount);

  // TODO: remove this unnecessary copy
  // Get shapesets
  const std::vector<const Shapeset *> &m_basesA =
      callVariant == TEST_TRIAL ? *m_testShapesets : *m_trialShapesets;
  const std::vector<const Shapeset *> &m_basesB =
      callVariant == TEST_TRIAL ? *m_trialShapesets : *m_testShapesets;
  std::vector<const Shapeset *> basesA(elementACount);
  for (int i = 0; i < elementACount; ++i)
    basesA[i] = m_basesA[elementIndicesA[i]];
  const Shapeset &basisB = *m_basesB[elementIndexB];

  // Find cached matrices; select integrators to calculate non-cached ones
  typedef std::pair<const Integrator *, const Shapeset *> QuadVariant;
  const QuadVariant CACHED(0, 0);
  std::vector<QuadVariant> quadVariants(elementACount);
  for (int i = 0; i < elementACount; ++i) {
    // Try to find matrix in cache
    const arma::Mat<ResultType> *cachedLocalWeakForm = 0;
    if (callVariant == TEST_TRIAL) {
      const int testElementIndex = elementIndicesA[i];
      const int trialElementIndex = elementIndexB;
      for (size_t n = 0; n < m_cache.extent(0); ++n)
        if (m_cache(n, trialElementIndex).first == testElementIndex) {
          cachedLocalWeakForm = &m_cache(n, trialElementIndex).second;
          break;
        }
    } else {
      const int testElementIndex = elementIndexB;
      const int trialElementIndex = elementIndicesA[i];
      for (size_t n = 0; n < m_cache.extent(0); ++n)
        if (m_cache(n, trialElementIndex).first == testElementIndex) {
          cachedLocalWeakForm = &m_cache(n, trialElementIndex).second;
          break;
        }
    }

    if (cachedLocalWeakForm) { // Matrix found in cache
      quadVariants[i] = CACHED;
      if (localDofIndexB == ALL_DOFS)
        result[i] = *cachedLocalWeakForm;
      else {
        if (callVariant == TEST_TRIAL)
          result[i] = cachedLocalWeakForm->col(localDofIndexB);
        else
          result[i] = cachedLocalWeakForm->row(localDofIndexB);
      }
    } else {
      const Integrator *integrator =
          callVariant == TEST_TRIAL
              ? &selectIntegrator(elementIndicesA[i], elementIndexB,
                                  nominalDistance)
              : &selectIntegrator(elementIndexB, elementIndicesA[i],
                                  nominalDistance);
      quadVariants[i] = QuadVariant(integrator, basesA[i]);
    }
  }

  // Integration will proceed in batches of test elements having the same
  // "quadrature variant", i.e. integrator and shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<int> activeElementIndicesA;
  activeElementIndicesA.reserve(elementACount);
  std::vector<arma::Mat<ResultType> *> activeLocalResults;
  activeLocalResults.reserve(elementACount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    if (activeQuadVariant == CACHED)
      continue;
    const Integrator &activeIntegrator = *it->first;
    const Shapeset &activeBasisA = *it->second;

    // Find all the test elements for which quadrature should proceed
    // according to the current quadrature variant
    activeElementIndicesA.clear();
    activeLocalResults.clear();
    for (int indexA = 0; indexA < elementACount; ++indexA)
      if (quadVariants[indexA] == activeQuadVariant) {
        activeElementIndicesA.push_back(elementIndicesA[indexA]);
        activeLocalResults.push_back(&result[indexA]);
      }

    // Integrate!
    activeIntegrator.integrate(callVariant, activeElementIndicesA,
                               elementIndexB, activeBasisA, basisB,
                               localDofIndexB, activeLocalResults);

    // // Distribute the just calculated integrals into the result array
    // // that will be returned to caller
    // int i = 0;
    // for (int indexA = 0; indexA < elementACount; ++indexA)
    //     if (quadVariants[indexA] == activeQuadVariant)
    //         result[indexA] = localResult.slice(i++);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    evaluateLocalWeakForms(const std::vector<int> &testElementIndices,
                           const std::vector<int> &trialElementIndices,
                           Fiber::_2dArray<arma::Mat<ResultType>> &result,
                           CoordinateType nominalDistance) {
  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;

  const int testElementCount = testElementIndices.size();
  const int trialElementCount = trialElementIndices.size();
  result.set_size(testElementCount, trialElementCount);

  // Find cached matrices; select integrators to calculate non-cached ones
  typedef boost::tuples::tuple<const Integrator *, const Shapeset *,
                               const Shapeset *> QuadVariant;
  const QuadVariant CACHED(0, 0, 0);
  Fiber::_2dArray<QuadVariant> quadVariants(testElementCount,
                                            trialElementCount);

  for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
    for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
      const int activeTestElementIndex = testElementIndices[testIndex];
      const int activeTrialElementIndex = trialElementIndices[trialIndex];
      // Try to find matrix in cache
      const arma::Mat<ResultType> *cachedLocalWeakForm = 0;
      for (size_t n = 0; n < m_cache.extent(0); ++n)
        if (m_cache(n, activeTrialElementIndex).first ==
            activeTestElementIndex) {
          cachedLocalWeakForm = &m_cache(n, activeTrialElementIndex).second;
          break;
        }

      if (cachedLocalWeakForm) { // Matrix found in cache
        quadVariants(testIndex, trialIndex) = CACHED;
        result(testIndex, trialIndex) = *cachedLocalWeakForm;
      } else {
        const Integrator *integrator =
            &selectIntegrator(activeTestElementIndex, activeTrialElementIndex,
                              nominalDistance);
        quadVariants(testIndex, trialIndex) =
            QuadVariant(integrator, (*m_testShapesets)[activeTestElementIndex],
                        (*m_trialShapesets)[activeTrialElementIndex]);
      }
    }

  // Integration will proceed in batches of element pairs having the same
  // "quadrature variant", i.e. integrator, test shapeset and trial shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<ElementIndexPair> activeElementPairs;
  std::vector<arma::Mat<ResultType> *> activeLocalResults;
  activeElementPairs.reserve(testElementCount * trialElementCount);
  activeLocalResults.reserve(testElementCount * trialElementCount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    if (activeQuadVariant == CACHED)
      continue;
    const Integrator &activeIntegrator = *it->template get<0>();
    const Shapeset &activeTestShapeset = *it->template get<1>();
    const Shapeset &activeTrialShapeset = *it->template get<2>();

    // Find all the element pairs for which quadrature should proceed
    // according to the current quadrature variant
    activeElementPairs.clear();
    activeLocalResults.clear();
    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
      for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
        if (quadVariants(testIndex, trialIndex) == activeQuadVariant) {
          activeElementPairs.push_back(ElementIndexPair(
              testElementIndices[testIndex], trialElementIndices[trialIndex]));
          activeLocalResults.push_back(&result(testIndex, trialIndex));
        }

    // Integrate!
    activeIntegrator.integrate(activeElementPairs, activeTestShapeset,
                               activeTrialShapeset, activeLocalResults);

    // // Distribute the just calculated integrals into the result array
    // // that will be returned to caller
    // int i = 0;
    // for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
    //     for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
    //         if (quadVariants(testIndex, trialIndex) == activeQuadVariant)
    //             result(testIndex, trialIndex) = localResult.slice(i++);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
typename DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::CoordinateType
DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::estimateRelativeScale(CoordinateType minDist) const {
  return m_kernels->estimateRelativeScale(minDist);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::cacheSingularLocalWeakForms() {
  ElementIndexPairSet elementIndexPairs;
  findPairsOfAdjacentElements(elementIndexPairs);
  cacheLocalWeakForms(elementIndexPairs);
}

/** \brief Fill \p pairs with the list of pairs of indices of elements
        sharing at least one vertex. */
template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::findPairsOfAdjacentElements(ElementIndexPairSet &pairs)
    const {
  pairs.clear();

  if (!testAndTrialGridsAreIdentical())
    return; // we assume that nonidentical grids are always disjoint

  const RawGridGeometry<CoordinateType> &rawGeometry = *m_testRawGeometry;

  const arma::Mat<CoordinateType> &vertices = rawGeometry.vertices();
  const arma::Mat<int> &elementCornerIndices =
      rawGeometry.elementCornerIndices();

  const int vertexCount = vertices.n_cols;
  const int elementCount = elementCornerIndices.n_cols;
  const int maxCornerCount = elementCornerIndices.n_rows;

  typedef std::vector<int> ElementIndexVector;
  // ith entry: set of elements sharing vertex number i
  std::vector<ElementIndexVector> elementsAdjacentToVertex(vertexCount);

  for (int e = 0; e < elementCount; ++e)
    for (int v = 0; v < maxCornerCount; ++v) {
      const int index = elementCornerIndices(v, e);
      if (index >= 0)
        elementsAdjacentToVertex[index].push_back(e);
    }

  // Loop over vertex indices
  for (int v = 0; v < vertexCount; ++v) {
    const ElementIndexVector &adjacentElements = elementsAdjacentToVertex[v];
    // Add to pairs each pair of elements adjacent to vertex v
    const int adjacentElementCount = adjacentElements.size();
    for (int e1 = 0; e1 < adjacentElementCount; ++e1)
      for (int e2 = 0; e2 < adjacentElementCount; ++e2)
        pairs.insert(
            ElementIndexPair(adjacentElements[e1], adjacentElements[e2]));
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::cacheLocalWeakForms(const ElementIndexPairSet &
                                              elementIndexPairs) {
  tbb::tick_count start = tbb::tick_count::now();

  if (elementIndexPairs.empty())
    return;

  if (m_verbosityLevel >= VerbosityLevel::DEFAULT)
    std::cout << "Precalculating singular integrals..." << std::endl;

  // Get the maximum number of neighbours a trial element can have.
  // This will be used to allocate a cache of correct size.
  // This loops assumes that elementIndexPairs are sorted after the trial
  // element index first.
  size_t maxNeighbourCount = 0;
  for (typename ElementIndexPairSet::const_iterator it =
           elementIndexPairs.begin();
       it != elementIndexPairs.end(); /* nothing */) {
    int curTrialElementIndex = it->second;
    size_t neighbourCount = 1;
    for (++it;
         it != elementIndexPairs.end() && it->second == curTrialElementIndex;
         ++it) {
      ++neighbourCount;
    }
    maxNeighbourCount = std::max(maxNeighbourCount, neighbourCount);
  }

  // Allocate cache and initialize all test element indices to INVALID_INDEX
  size_t trialElementCount = m_trialRawGeometry->elementCount();
  m_cache.set_size(maxNeighbourCount, trialElementCount);
  for (size_t trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
    for (size_t nborIndex = 0; nborIndex < maxNeighbourCount; ++nborIndex)
      m_cache(nborIndex, trialIndex).first = INVALID_INDEX;

  // Find cached matrices; select integrators to calculate non-cached ones
  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;
  typedef boost::tuples::tuple<const Integrator *, const Shapeset *,
                               const Shapeset *> QuadVariant;
  const int elementPairCount = elementIndexPairs.size();
  std::vector<QuadVariant> quadVariants(elementPairCount);

  typedef typename ElementIndexPairSet::const_iterator ElementIndexPairIterator;
  typedef typename std::vector<QuadVariant>::iterator QuadVariantIterator;
  {
    ElementIndexPairIterator pairIt = elementIndexPairs.begin();
    QuadVariantIterator qvIt = quadVariants.begin();
    for (; pairIt != elementIndexPairs.end(); ++pairIt, ++qvIt) {
      const int testElementIndex = pairIt->first;
      const int trialElementIndex = pairIt->second;
      const Integrator *integrator =
          &selectIntegrator(testElementIndex, trialElementIndex);
      *qvIt = QuadVariant(integrator, (*m_testShapesets)[testElementIndex],
                          (*m_trialShapesets)[trialElementIndex]);
    }
  }

  // Integration will proceed in batches of element pairs having the same
  // "quadrature variant", i.e. integrator, test shapeset and trial shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<ElementIndexPair> activeElementPairs;
  std::vector<arma::Mat<ResultType> *> activeLocalResults;
  activeElementPairs.reserve(elementPairCount);
  activeLocalResults.reserve(elementPairCount);
  // m_cache.rehash(int(elementIndexPairs.size() / m_cache.max_load_factor() +
  // 1));

  int maxThreadCount = 1;
  if (!m_parallelizationOptions.isOpenClEnabled()) {
    if (m_parallelizationOptions.maxThreadCount() ==
        ParallelizationOptions::AUTO)
      maxThreadCount = tbb::task_scheduler_init::automatic;
    else
      maxThreadCount = m_parallelizationOptions.maxThreadCount();
  }
  tbb::task_scheduler_init scheduler(maxThreadCount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    const Integrator &activeIntegrator = *it->template get<0>();
    const Shapeset &activeTestShapeset = *it->template get<1>();
    const Shapeset &activeTrialShapeset = *it->template get<2>();

    // Find all the element pairs for which quadrature should proceed
    // according to the current quadrature variant
    activeElementPairs.clear();
    activeLocalResults.clear();
    {
      ElementIndexPairIterator pairIt = elementIndexPairs.begin();
      QuadVariantIterator qvIt = quadVariants.begin();
      size_t neighbourIndex = 0;                    // Cache row index
      size_t curTrialElementIndex = pairIt->second; // Cache column index
      for (; pairIt != elementIndexPairs.end(); ++pairIt, ++qvIt) {
        if (pairIt->second != curTrialElementIndex) {
          // Go to the next column of cache
          curTrialElementIndex = pairIt->second;
          neighbourIndex = 0;
        }
        if (*qvIt == activeQuadVariant) {
          activeElementPairs.push_back(*pairIt);
          m_cache(neighbourIndex, curTrialElementIndex).first = pairIt->first;
          activeLocalResults.push_back(
              &m_cache(neighbourIndex, curTrialElementIndex).second);
        }
        // Increment cache row index
        ++neighbourIndex;
      }
    }

    // Integrate!
    // Old serial version
    // activeIntegrator.integrate(activeElementPairs, activeTestShapeset,
    //                            activeTrialShapeset, localResult);

    typedef SingularIntegralCalculatorLoopBody<BasisFunctionType, KernelType,
                                               ResultType> Body;
    {
      Fiber::SerialBlasRegion region;
      tbb::parallel_for(
          tbb::blocked_range<size_t>(0, activeElementPairs.size()),
          Body(activeIntegrator, activeElementPairs, activeTestShapeset,
               activeTrialShapeset, activeLocalResults));
    }
  }
  tbb::tick_count end = tbb::tick_count::now();
  if (m_verbosityLevel >= VerbosityLevel::DEFAULT)
    std::cout << "Precalculation of singular integrals took "
              << (end - start).seconds() << " s" << std::endl;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
const TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> &
DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::selectIntegrator(int testElementIndex,
                                       int trialElementIndex,
                                       CoordinateType nominalDistance) {
  DoubleQuadratureDescriptor desc = m_quadDescSelector->quadratureDescriptor(
      testElementIndex, trialElementIndex, nominalDistance);
  return getIntegrator(desc);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
const TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> &
DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::getIntegrator(const DoubleQuadratureDescriptor &desc) {
  typename IntegratorMap::iterator it = m_testKernelTrialIntegrators.find(desc);
  // Note: as far as I understand TBB's docs, .end() keeps pointing to the
  // same element even if another thread inserts a new element into the map
  if (it == m_testKernelTrialIntegrators.end()) {
    tbb::mutex::scoped_lock lock(m_integratorCreationMutex);
    it = m_testKernelTrialIntegrators.find(desc);
    if (it == m_testKernelTrialIntegrators.end()) {
      // std::cout << "getIntegrator(" << desc
      //           << "): creating an integrator" << std::endl;
      // Integrator doesn't exist yet and must be created.
      arma::Mat<CoordinateType> testPoints, trialPoints;
      std::vector<CoordinateType> testWeights, trialWeights;
      bool isTensor;
      m_quadRuleFamily->fillQuadraturePointsAndWeights(
          desc, testPoints, trialPoints, testWeights, trialWeights, isTensor);
      Integrator *integrator = 0;
      if (isTensor) {
        typedef SeparableNumericalTestKernelTrialIntegrator<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
        ConcreteIntegrator;
        integrator = new ConcreteIntegrator(
            testPoints, trialPoints, testWeights, trialWeights,
            *m_testGeometryFactory, *m_trialGeometryFactory, *m_testRawGeometry,
            *m_trialRawGeometry, *m_testTransformations, *m_kernels,
            *m_trialTransformations, *m_integral, *m_openClHandler);
      } else {
        typedef NonseparableNumericalTestKernelTrialIntegrator<
            BasisFunctionType, KernelType, ResultType, GeometryFactory>
        ConcreteIntegrator;
        integrator = new ConcreteIntegrator(
            testPoints, trialPoints, testWeights, *m_testGeometryFactory,
            *m_trialGeometryFactory, *m_testRawGeometry, *m_trialRawGeometry,
            *m_testTransformations, *m_kernels, *m_trialTransformations,
            *m_integral, *m_openClHandler);
      }

      // Attempt to insert the newly created integrator into the map
      std::pair<typename IntegratorMap::iterator, bool> result =
          m_testKernelTrialIntegrators.insert(std::make_pair(desc, integrator));
      if (result.second)
        // Insertion succeeded. The newly created integrator will be deleted in
        // our own destructor
        ;
      else
        // Insertion failed -- another thread was faster. Delete the newly
        // created integrator.
        delete integrator;

      //    if (result.second)
      //        std::cout << "getIntegrator(: " << desc << "): insertion
      // succeeded" << std::endl;
      //    else
      //        std::cout << "getIntegrator(: " << desc << "): insertion failed"
      // << std::endl;

      // Return pointer to the integrator that ended up in the map.
      it = result.first;
    }
  }
  return *it->second;
}

} // namespace Fiber
