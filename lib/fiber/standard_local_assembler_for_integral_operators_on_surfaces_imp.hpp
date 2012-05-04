// Copyright (C) 2011-2012 by the Fiber Authors
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

// Keep IDEs happy
#include "standard_local_assembler_for_integral_operators_on_surfaces.hpp"

#include "nonseparable_numerical_test_kernel_trial_integrator.hpp"
#include "separable_numerical_test_kernel_trial_integrator.hpp"

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace Fiber
{

namespace
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class SingularIntegralCalculatorLoopBody
{
public:
    typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> Integrator;
    typedef typename Integrator::ElementIndexPair ElementIndexPair;

    SingularIntegralCalculatorLoopBody(
            const Integrator& activeIntegrator,
            const std::vector<ElementIndexPair>& activeElementPairs,
            const Basis<BasisFunctionType>& activeTestBasis,
            const Basis<BasisFunctionType>& activeTrialBasis,
            arma::Cube<ResultType>& localResult) :
        m_activeIntegrator(activeIntegrator),
        m_activeElementPairs(activeElementPairs),
        m_activeTestBasis(activeTestBasis),
        m_activeTrialBasis(activeTrialBasis),
        m_localResult(localResult) {
    }

    void operator() (const tbb::blocked_range<size_t>& r) const {
        // copy the relevant subset of m_activeElementPairs into
        // localActiveElementPairs
        std::vector<ElementIndexPair> localActiveElementPairs(
                    &m_activeElementPairs[r.begin()],
                    &m_activeElementPairs[r.end()]);
        arma::Cube<ResultType> localLocalResult(
                    &m_localResult(0, 0, r.begin()),
                    m_localResult.n_rows,
                    m_localResult.n_cols,
                    r.size(),
                    false /* copy_aux_mem */);
        m_activeIntegrator.integrate(localActiveElementPairs, m_activeTestBasis,
                                     m_activeTrialBasis, localLocalResult);
    }

private:
    const Integrator& m_activeIntegrator;
    const std::vector<ElementIndexPair>& m_activeElementPairs;
    const Basis<BasisFunctionType>& m_activeTestBasis;
    const Basis<BasisFunctionType>& m_activeTrialBasis;
    mutable arma::Cube<ResultType>& m_localResult;
};

} // namespace

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
StandardLocalAssemblerForIntegralOperatorsOnSurfaces(
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const std::vector<const Basis<BasisFunctionType>*>& testBases,
        const std::vector<const Basis<BasisFunctionType>*>& trialBases,
        const Expression<CoordinateType>& testExpression,
        const Kernel<KernelType>& kernel,
        const Expression<CoordinateType>& trialExpression,
        const OpenClHandler& openClHandler,
        const ParallelisationOptions& parallelisationOptions,
        bool cacheSingularIntegrals,
        const AccuracyOptions& accuracyOptions) :
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testBases(testBases),
    m_trialBases(trialBases),
    m_testExpression(testExpression),
    m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_openClHandler(openClHandler),
    m_parallelisationOptions(parallelisationOptions),
    m_accuracyOptions(accuracyOptions)
{
    if (rawGeometry.vertices().n_rows != 3)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces::"
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces(): "
                "vertex coordinates must be three-dimensional");
    if (rawGeometry.elementCornerIndices().n_rows < 3 ||
            4 < rawGeometry.elementCornerIndices().n_rows)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces::"
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces(): "
                "all elements must be triangular or quadrilateral");
    const int elementCount = rawGeometry.elementCornerIndices().n_cols;
    if (!rawGeometry.auxData().is_empty() &&
            rawGeometry.auxData().n_cols != elementCount)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces::"
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces(): "
                "number of columns of auxData must match that of "
                "elementCornerIndices");
    if (testBases.size() != elementCount)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces::"
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces(): "
                "size of testBases must match the number of columns of "
                "elementCornerIndices");
    if (trialBases.size() != elementCount)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces::"
                "StandardLocalAssemblerForIntegralOperatorsOnSurfaces(): "
                "size of trialBases must match the number of columns of "
                "elementCornerIndices");

    if (cacheSingularIntegrals)
        cacheSingularLocalWeakForms();
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
~StandardLocalAssemblerForIntegralOperatorsOnSurfaces()
{
    // Note: obviously the destructor is assumed to be called only after
    // all threads have ceased using the assembler!

    for (typename IntegratorMap::const_iterator it = m_TestKernelTrialIntegrators.begin();
         it != m_TestKernelTrialIntegrators.end(); ++it)
        delete it->second;
    m_TestKernelTrialIntegrators.clear();
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
evaluateLocalWeakForms(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        LocalDofIndex localDofIndexB,
        std::vector<arma::Mat<ResultType> >& result)
{
    typedef Basis<BasisFunctionType> Basis;

    const int elementACount = elementIndicesA.size();
    result.resize(elementACount);

    // TODO: remove this unnecessary copy
    // Get bases
    const std::vector<const Basis*>& m_basesA =
            callVariant == TEST_TRIAL ? m_testBases : m_trialBases;
    const std::vector<const Basis*>& m_basesB =
            callVariant == TEST_TRIAL ? m_trialBases : m_testBases;
    std::vector<const Basis*> basesA(elementACount);
    for (int i = 0; i < elementACount; ++i)
        basesA[i] = m_basesA[elementIndicesA[i]];
    const Basis& basisB = *m_basesB[elementIndexB];

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    const QuadVariant CACHED(0, 0);
    std::vector<QuadVariant> quadVariants(elementACount);
    for (int i = 0; i < elementACount; ++i) {
        typename Cache::const_iterator it = m_cache.find(
                    callVariant == TEST_TRIAL ?
                        ElementIndexPair(elementIndicesA[i], elementIndexB) :
                        ElementIndexPair(elementIndexB, elementIndicesA[i]));
        if (it != m_cache.end()) { // Matrix found in cache
            quadVariants[i] = CACHED;
            if (localDofIndexB == ALL_DOFS)
                result[i] = it->second;
            else {
                if (callVariant == TEST_TRIAL)
                    result[i] = it->second.col(localDofIndexB);
                else
                    result[i] = it->second.row(localDofIndexB);
            }
        } else {
            const Integrator* integrator =
                    callVariant == TEST_TRIAL ?
                        &selectIntegrator(elementIndicesA[i], elementIndexB) :
                        &selectIntegrator(elementIndexB, elementIndicesA[i]);
            quadVariants[i] = QuadVariant(integrator, basesA[i]);
        }
    }

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and basis

    // Find all the unique quadrature variants present
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<int> activeElementIndicesA;
    activeElementIndicesA.reserve(elementACount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it) {
        const QuadVariant activeQuadVariant = *it;
        if (activeQuadVariant == CACHED)
            continue;
        const Integrator& activeIntegrator = *it->first;
        const Basis& activeBasisA = *it->second;

        // Find all the test elements for which quadrature should proceed
        // according to the current quadrature variant
        activeElementIndicesA.clear();
        for (int indexA = 0; indexA < elementACount; ++indexA)
            if (quadVariants[indexA] == activeQuadVariant)
                activeElementIndicesA.push_back(elementIndicesA[indexA]);

        // Integrate!
        arma::Cube<ResultType> localResult;
        activeIntegrator.integrate(callVariant,
                                   activeElementIndicesA, elementIndexB,
                                   activeBasisA, basisB, localDofIndexB,
                                   localResult);

        // Distribute the just calculated integrals into the result array
        // that will be returned to caller
        int i = 0;
        for (int indexA = 0; indexA < elementACount; ++indexA)
            if (quadVariants[indexA] == activeQuadVariant)
                result[indexA] = localResult.slice(i++);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& testElementIndices,
        const std::vector<int>& trialElementIndices,
        Fiber::Array2d<arma::Mat<ResultType> >& result)
{
    typedef Fiber::Basis<BasisFunctionType> Basis;

    const int testElementCount = testElementIndices.size();
    const int trialElementCount = trialElementIndices.size();
    result.set_size(testElementCount, trialElementCount);

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef boost::tuples::tuple<const Integrator*, const Basis*, const Basis*>
            QuadVariant;
    const QuadVariant CACHED(0, 0, 0);
    Fiber::Array2d<QuadVariant> quadVariants(testElementCount, trialElementCount);

    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
            const int activeTestElementIndex = testElementIndices[testIndex];
            const int activeTrialElementIndex = trialElementIndices[trialIndex];
            typename Cache::const_iterator it = m_cache.find(
                        ElementIndexPair(activeTestElementIndex,
                                         activeTrialElementIndex));
            if (it != m_cache.end()) { // Matrix found in cache
                quadVariants(testIndex, trialIndex) = CACHED;
                result(testIndex, trialIndex) = it->second;
            } else {
                const Integrator* integrator =
                        &selectIntegrator(activeTestElementIndex,
                                          activeTrialElementIndex);
                quadVariants(testIndex, trialIndex) = QuadVariant(
                            integrator, m_testBases[activeTestElementIndex],
                            m_trialBases[activeTrialElementIndex]);
            }
        }

    // Integration will proceed in batches of element pairs having the same
    // "quadrature variant", i.e. integrator, test basis and trial basis

    // Find all the unique quadrature variants present
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<ElementIndexPair> activeElementPairs;
    activeElementPairs.reserve(testElementCount * trialElementCount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it) {
        const QuadVariant activeQuadVariant = *it;
        if (activeQuadVariant == CACHED)
            continue;
        const Integrator& activeIntegrator = *it->template get<0>();
        const Basis& activeTestBasis  = *it->template get<1>();
        const Basis& activeTrialBasis = *it->template get<2>();

        // Find all the element pairs for which quadrature should proceed
        // according to the current quadrature variant
        activeElementPairs.clear();
        for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
            for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
                if (quadVariants(testIndex, trialIndex) == activeQuadVariant)
                    activeElementPairs.push_back(
                                ElementIndexPair(testElementIndices[testIndex],
                                                 trialElementIndices[trialIndex]));

        // Integrate!
        arma::Cube<ResultType> localResult;
        activeIntegrator.integrate(activeElementPairs, activeTestBasis,
                                   activeTrialBasis, localResult);

        // Distribute the just calculated integrals into the result array
        // that will be returned to caller
        int i = 0;
        for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
            for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
                if (quadVariants(testIndex, trialIndex) == activeQuadVariant)
                    result(testIndex, trialIndex) = localResult.slice(i++);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& elementIndices,
        std::vector<arma::Mat<ResultType> >& result)
{
    // This overload is mostly useful only for the identity operator
    throw std::runtime_error("StandardLocalAssemblerForIntegralOperatorsOnSurfaces::"
                             "evaluateLocalWeakForms(): "
                             "this overload not implemented yet");
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
cacheSingularLocalWeakForms()
{
    ElementIndexPairSet elementIndexPairs;
    findPairsOfAdjacentElements(elementIndexPairs);
    cacheLocalWeakForms(elementIndexPairs);
}

/** \brief Fill \p pairs with the list of pairs of indices of elements
        sharing at least one vertex. */
template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
findPairsOfAdjacentElements(ElementIndexPairSet& pairs) const
{
    const arma::Mat<CoordinateType>& vertices = m_rawGeometry.vertices();
    const arma::Mat<int>& elementCornerIndices =
            m_rawGeometry.elementCornerIndices();

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

    pairs.clear();
    // Loop over vertex indices
    for (int v = 0; v < vertexCount; ++v) {
        const ElementIndexVector& adjacentElements = elementsAdjacentToVertex[v];
        // Add to pairs each pair of elements adjacent to vertex v
        const int adjacentElementCount = adjacentElements.size();
        for (int e1 = 0; e1 < adjacentElementCount; ++e1)
            for (int e2 = 0; e2 < adjacentElementCount; ++e2)
                pairs.insert(ElementIndexPair(adjacentElements[e1],
                                              adjacentElements[e2]));
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
cacheLocalWeakForms(const ElementIndexPairSet& elementIndexPairs)
{
    typedef Fiber::Basis<BasisFunctionType> Basis;

    const int elementPairCount = elementIndexPairs.size();

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef boost::tuples::tuple<const Integrator*, const Basis*, const Basis*>
            QuadVariant;
    std::vector<QuadVariant> quadVariants(elementPairCount);

    typedef typename ElementIndexPairSet::const_iterator
            ElementIndexPairIterator;
    typedef typename std::vector<QuadVariant>::iterator QuadVariantIterator;
    {
        ElementIndexPairIterator pairIt = elementIndexPairs.begin();
        QuadVariantIterator qvIt = quadVariants.begin();
        for (; pairIt != elementIndexPairs.end(); ++pairIt, ++qvIt) {
            const int testElementIndex = pairIt->first;
            const int trialElementIndex = pairIt->second;
            const Integrator* integrator =
                    &selectIntegrator(testElementIndex, trialElementIndex);
            *qvIt = QuadVariant(integrator,
                                m_testBases[testElementIndex],
                                m_trialBases[trialElementIndex]);
        }
    }

    // Integration will proceed in batches of element pairs having the same
    // "quadrature variant", i.e. integrator, test basis and trial basis

    // Find all the unique quadrature variants present
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<ElementIndexPair> activeElementPairs;
    activeElementPairs.reserve(elementPairCount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it) {
        const QuadVariant activeQuadVariant = *it;
        const Integrator& activeIntegrator = *it->template get<0>();
        const Basis& activeTestBasis  = *it->template get<1>();
        const Basis& activeTrialBasis = *it->template get<2>();

        // Find all the element pairs for which quadrature should proceed
        // according to the current quadrature variant
        activeElementPairs.clear();
        {
            ElementIndexPairIterator pairIt = elementIndexPairs.begin();
            QuadVariantIterator qvIt = quadVariants.begin();
            for (; pairIt != elementIndexPairs.end(); ++pairIt, ++qvIt)
                if (*qvIt == activeQuadVariant)
                    activeElementPairs.push_back(*pairIt);
        }

        // Integrate!
        arma::Cube<ResultType> localResult(activeTestBasis.size(),
                                          activeTrialBasis.size(),
                                          activeElementPairs.size());
        // Old serial version
        // activeIntegrator.integrate(activeElementPairs, activeTestBasis,
        //                            activeTrialBasis, localResult);

        const int maxThreadCount =
                m_parallelisationOptions.mode() == ParallelisationOptions::TBB ?
                    m_parallelisationOptions.maxThreadCount() : 1;
        tbb::task_scheduler_init scheduler(maxThreadCount);
        typedef SingularIntegralCalculatorLoopBody<
                BasisFunctionType, KernelType, ResultType> Body;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, activeElementPairs.size()),
                          Body(activeIntegrator,
                               activeElementPairs, activeTestBasis, activeTrialBasis,
                               localResult));

        {
            ElementIndexPairIterator pairIt = elementIndexPairs.begin();
            QuadVariantIterator qvIt = quadVariants.begin();
            int i = 0;
            for (; pairIt != elementIndexPairs.end(); ++pairIt, ++qvIt)
                if (*qvIt == activeQuadVariant)
                    m_cache[*pairIt] = localResult.slice(i++);
        }
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
const TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>&
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
selectIntegrator(int testElementIndex, int trialElementIndex)
{
    DoubleQuadratureDescriptor desc;

    // Get corner indices of the specified elements
    arma::Col<int> testElementCornerIndices =
            m_rawGeometry.elementCornerIndices(testElementIndex);
    arma::Col<int> trialElementCornerIndices =
            m_rawGeometry.elementCornerIndices(trialElementIndex);

    desc.topology = determineElementPairTopologyIn3D(
                testElementCornerIndices, trialElementCornerIndices);

    if (desc.topology.type == ElementPairTopology::Disjoint) {
        desc.testOrder = regularOrder(testElementIndex, TEST);
        desc.trialOrder = regularOrder(testElementIndex, TRIAL);
    } else { // singular integral
        desc.testOrder = singularOrder(testElementIndex, TEST);
        desc.trialOrder = singularOrder(testElementIndex, TRIAL);
    }

    return getIntegrator(desc);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
regularOrder(int elementIndex, ElementType elementType) const
{
    // TODO:
    // 1. Check the size of elements and the distance between them
    //    and estimate the variability of the kernel
    // 2. Take into account the fact that elements might be isoparametric.

    const QuadratureOptions& options = m_accuracyOptions.doubleRegular;

    if (options.mode == QuadratureOptions::EXACT_ORDER)
        return options.order;
    else {
        // Order required for exact quadrature on affine elements with a constant kernel
        int elementOrder = (elementType == TEST ?
                                m_testBases[elementIndex]->order() :
                                m_trialBases[elementIndex]->order());
        int minimumOrder = ((elementOrder + 1) + 1) / 2;
        return minimumOrder + options.orderIncrement;
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
singularOrder(int elementIndex, ElementType elementType) const
{
    // TODO:
    // 1. Check the size of elements and estimate the variability of the
    //    (Sauter-Schwab-transformed) kernel
    // 2. Take into account the fact that elements might be isoparametric.

    const QuadratureOptions& options = m_accuracyOptions.doubleSingular;

    if (options.mode == QuadratureOptions::EXACT_ORDER)
        return options.order;
    else {
        // Order required for exact quadrature on affine elements with a constant kernel
        int elementOrder = (elementType == TEST ?
                                m_testBases[elementIndex]->order() :
                                m_trialBases[elementIndex]->order());
        int minimumOrder = ((elementOrder + 1) + 1) / 2;
        return minimumOrder + 3 + options.orderIncrement;
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
const TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>&
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
getIntegrator(const DoubleQuadratureDescriptor& desc)
{
    typename IntegratorMap::const_iterator it = m_TestKernelTrialIntegrators.find(desc);
    // Note: as far as I understand TBB's docs, .end() keeps pointing to the
    // same element even if another thread inserts a new element into the map
    if (it != m_TestKernelTrialIntegrators.end()) {
        //std::cout << "getIntegrator(: " << desc << "): integrator found" << std::endl;
        return *it->second;
    }
    //std::cout << "getIntegrator(: " << desc << "): integrator not found" << std::endl;

    // Integrator doesn't exist yet and must be created.
    Integrator* integrator = 0;
    const ElementPairTopology& topology = desc.topology;
    if (topology.type == ElementPairTopology::Disjoint) {
        // Create a tensor rule
        arma::Mat<CoordinateType> testPoints, trialPoints;
        std::vector<CoordinateType> testWeights, trialWeights;

        fillSingleQuadraturePointsAndWeights(topology.testVertexCount,
                                             desc.testOrder,
                                             testPoints, testWeights);
        fillSingleQuadraturePointsAndWeights(topology.trialVertexCount,
                                             desc.trialOrder,
                                             trialPoints, trialWeights);
        typedef SeparableNumericalTestKernelTrialIntegrator<BasisFunctionType,
                KernelType, ResultType, GeometryFactory> ConcreteIntegrator;
        integrator = new ConcreteIntegrator(
                    testPoints, trialPoints, testWeights, trialWeights,
                    m_geometryFactory, m_rawGeometry,
                    m_testExpression, m_kernel, m_trialExpression,
                    m_openClHandler);
    } else {
        arma::Mat<CoordinateType> testPoints, trialPoints;
        std::vector<CoordinateType> weights;

        fillDoubleSingularQuadraturePointsAndWeights(
                    desc, testPoints, trialPoints, weights);
        typedef NonseparableNumericalTestKernelTrialIntegrator<BasisFunctionType,
                KernelType, ResultType, GeometryFactory> ConcreteIntegrator;
        integrator = new ConcreteIntegrator(
                    testPoints, trialPoints, weights,
                    m_geometryFactory, m_rawGeometry,
                    m_testExpression, m_kernel, m_trialExpression,
                    m_openClHandler);
    }

    // Attempt to insert the newly created integrator into the map
    std::pair<typename IntegratorMap::iterator, bool> result =
            m_TestKernelTrialIntegrators.insert(std::make_pair(desc, integrator));
    if (result.second)
        // Insertion succeeded. The newly created integrator will be deleted in
        // our own destructor
        ;
    else
        // Insertion failed -- another thread was faster. Delete the newly
        // created integrator.
        delete integrator;

    //    if (result.second)
    //        std::cout << "getIntegrator(: " << desc << "): insertion succeeded" << std::endl;
    //    else
    //        std::cout << "getIntegrator(: " << desc << "): insertion failed" << std::endl;

    // Return pointer to the integrator that ended up in the map.
    return *result.first->second;
}

} // namespace Fiber
