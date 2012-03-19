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

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
StandardLocalAssemblerForIntegralOperatorsOnSurfaces(
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<ValueType>& rawGeometry,
        const std::vector<const Basis<ValueType>*>& testBases,
        const std::vector<const Basis<ValueType>*>& trialBases,
        const Expression<ValueType>& testExpression,
        const Kernel<ValueType>& kernel,
        const Expression<ValueType>& trialExpression,
        const OpenClHandler& openClHandler,
        bool cacheSingularIntegrals) :
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testBases(testBases),
    m_trialBases(trialBases),
    m_testExpression(testExpression),
    m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_openClHandler(openClHandler)
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

template <typename ValueType, typename GeometryFactory>
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
~StandardLocalAssemblerForIntegralOperatorsOnSurfaces()
{
    // Note: obviously the destructor is assumed to be called only after
    // all threads have ceased using the assembler!

    for (typename IntegratorMap::const_iterator it = m_TestKernelTrialIntegrators.begin();
         it != m_TestKernelTrialIntegrators.end(); ++it)
        delete it->second;
    m_TestKernelTrialIntegrators.clear();
}

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
evaluateLocalWeakForms(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        LocalDofIndex localDofIndexB,
        std::vector<arma::Mat<ValueType> >& result)
{
    typedef TestKernelTrialIntegrator<ValueType> Integrator;
    typedef Basis<ValueType> Basis;

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

    // TODO: possibly precalculate diagonal elements

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    const QuadVariant CACHED(0, 0);
    std::vector<QuadVariant> quadVariants(elementACount);
    for (int i = 0; i < elementACount; ++i)
    {
        typename Cache::const_iterator it = m_cache.find(
                    ElementIndexPair(elementIndicesA[i], elementIndexB));
        if (it != m_cache.end()) // Matrix found in cache
        {
            quadVariants[i] = CACHED;
            if (localDofIndexB == ALL_DOFS)
                result[i] = it->second;
            else
            {
                result[i].set_size(it->second.n_rows, 1);
                result[i] = it->second.col(localDofIndexB);
            }
        }
        else
        {
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
         it != uniqueQuadVariants.end(); ++it)
    {
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
        arma::Cube<ValueType> localResult;
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

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& testElementIndices,
        const std::vector<int>& trialElementIndices,
        Fiber::Array2D<arma::Mat<ValueType> >& result)
{
    typedef Fiber::TestKernelTrialIntegrator<ValueType> Integrator;
    typedef Fiber::Basis<ValueType> Basis;

    const int testElementCount = testElementIndices.size();
    const int trialElementCount = trialElementIndices.size();
    result.set_size(testElementCount, trialElementCount);

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef boost::tuples::tuple<const Integrator*, const Basis*, const Basis*>
            QuadVariant;
    const QuadVariant CACHED(0, 0, 0);
    Fiber::Array2D<QuadVariant> quadVariants(testElementCount, trialElementCount);

    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
        {
            const int activeTestElementIndex = testElementIndices[testIndex];
            const int activeTrialElementIndex = trialElementIndices[trialIndex];
            typename Cache::const_iterator it = m_cache.find(
                    ElementIndexPair(activeTestElementIndex,
                                     activeTrialElementIndex));
            if (it != m_cache.end()) // Matrix found in cache
            {
                quadVariants(testIndex, trialIndex) = CACHED;
                result(testIndex, trialIndex) = it->second;
            }
            else
            {
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
         it != uniqueQuadVariants.end(); ++it)
    {
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
        arma::Cube<ValueType> localResult;
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

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
cacheSingularLocalWeakForms()
{
    ElementIndexPairSet elementIndexPairs;
    findPairsOfAdjacentElements(elementIndexPairs);
    cacheLocalWeakForms(elementIndexPairs);
}

/** \brief Fill \p pairs with the list of pairs of indices of elements
        sharing at least one vertex. */
template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
findPairsOfAdjacentElements(ElementIndexPairSet& pairs) const
{
    const arma::Mat<ValueType>& vertices = m_rawGeometry.vertices();
    const arma::Mat<int>& elementCornerIndices =
            m_rawGeometry.elementCornerIndices();

    const int vertexCount = vertices.n_cols;
    const int elementCount = elementCornerIndices.n_cols;
    const int maxCornerCount = elementCornerIndices.n_rows;

    typedef std::vector<int> ElementIndexVector;
    // ith entry: set of elements sharing vertex number i
    std::vector<ElementIndexVector> elementsAdjacentToVertex(vertexCount);

    for (int e = 0; e < elementCount; ++e)
        for (int v = 0; v < maxCornerCount; ++v)
        {
            const int index = elementCornerIndices(v, e);
            if (index >= 0)
                elementsAdjacentToVertex[index].push_back(e);
        }

    pairs.clear();
    // Loop over vertex indices
    for (int v = 0; v < vertexCount; ++v)
    {
        const ElementIndexVector& adjacentElements = elementsAdjacentToVertex[v];
        // Add to pairs each pair of elements adjacent to vertex v
        const int adjacentElementCount = adjacentElements.size();
        for (int e1 = 0; e1 < adjacentElementCount; ++e1)
            for (int e2 = e1; e2 < adjacentElementCount; ++e2)
                pairs.insert(ElementIndexPair(adjacentElements[e1],
                                              adjacentElements[e2]));
    }
}

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
cacheLocalWeakForms(const ElementIndexPairSet& elementIndexPairs)
{
    typedef Fiber::TestKernelTrialIntegrator<ValueType> Integrator;
    typedef Fiber::Basis<ValueType> Basis;

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
        for (; pairIt != elementIndexPairs.end(); ++pairIt, ++qvIt)
        {
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
         it != uniqueQuadVariants.end(); ++it)
    {
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
        arma::Cube<ValueType> localResult;
        activeIntegrator.integrate(activeElementPairs, activeTestBasis,
                                   activeTrialBasis, localResult);

        // Distribute the just calculated integrals into the result array
        // that will be returned to caller
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

template <typename ValueType, typename GeometryFactory>
const TestKernelTrialIntegrator<ValueType>&
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
selectIntegrator(int testElementIndex, int trialElementIndex) {
    DoubleQuadratureDescriptor desc;

    // Get corner indices of the specified elements
    arma::Col<int> testElementCornerIndices =
            m_rawGeometry.elementCornerIndices(testElementIndex);
    arma::Col<int> trialElementCornerIndices =
            m_rawGeometry.elementCornerIndices(trialElementIndex);

    desc.topology = determineElementPairTopologyIn3D(
                testElementCornerIndices, trialElementCornerIndices);

    desc.testOrder = m_testBases[testElementIndex]->order();
    desc.trialOrder = m_trialBases[trialElementIndex]->order();

    if (desc.topology.type == ElementPairTopology::Disjoint)
    {
        desc.testOrder +=
                regularOrderIncrement(testElementIndex);
        desc.trialOrder +=
                regularOrderIncrement(trialElementIndex);
    }
    else // singular integral
    {
        desc.testOrder +=
                singularOrderIncrement(testElementIndex);
        desc.trialOrder +=
                singularOrderIncrement(trialElementIndex);
    }

    return getIntegrator(desc);
}

template <typename ValueType, typename GeometryFactory>
int
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
regularOrderIncrement(int elementIndex) const
{
    // Note: this function will make use of options supplied to the integrator
    // in its constructor

    // TODO:
    // 1. Check the size of elements and the distance between them
    //    and estimate the variability of the kernel
    // 2. Take into account the fact that elements might be isoparametric.

    // Make quadrature exact for a constant kernel and affine elements
    return 4;
}

template <typename ValueType, typename GeometryFactory>
int
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
singularOrderIncrement(int elementIndex) const
{
    // Note: this function will make use of options supplied to the integrator
    // in its constructor

    // TODO:
    // 1. Check the size of elements and estimate the variability of the
    //    (Sauter-Schwab-transformed) kernel
    // 2. Take into account the fact that elements might be isoparametric.

    // Make quadrature exact for a constant kernel and affine elements
    return 4;
}

template <typename ValueType, typename GeometryFactory>
const TestKernelTrialIntegrator<ValueType>&
StandardLocalAssemblerForIntegralOperatorsOnSurfaces<ValueType, GeometryFactory>::
getIntegrator(const DoubleQuadratureDescriptor& desc)
{
    typename IntegratorMap::const_iterator it = m_TestKernelTrialIntegrators.find(desc);
    // Note: as far as I understand TBB's docs, .end() keeps pointing to the
    // same element even if another thread inserts a new element into the map
    if (it != m_TestKernelTrialIntegrators.end())
    {
        //std::cout << "getIntegrator(: " << desc << "): integrator found" << std::endl;
        return *it->second;
    }
    //std::cout << "getIntegrator(: " << desc << "): integrator not found" << std::endl;

    // Integrator doesn't exist yet and must be created.
    TestKernelTrialIntegrator<ValueType>* integrator = 0;
    const ElementPairTopology& topology = desc.topology;
    if (topology.type == ElementPairTopology::Disjoint)
    {
        // Create a tensor rule
        arma::Mat<ValueType> testPoints, trialPoints;
        std::vector<ValueType> testWeights, trialWeights;

        fillSingleQuadraturePointsAndWeights(topology.testVertexCount,
                                             desc.testOrder,
                                             testPoints, testWeights);
        fillSingleQuadraturePointsAndWeights(topology.trialVertexCount,
                                             desc.trialOrder,
                                             trialPoints, trialWeights);
        typedef SeparableNumericalTestKernelTrialIntegrator<ValueType, GeometryFactory> Integrator;
        integrator = new Integrator(
                        testPoints, trialPoints, testWeights, trialWeights,
                        m_geometryFactory, m_rawGeometry,
                        m_testExpression, m_kernel, m_trialExpression,
                        m_openClHandler);
    }
    else
    {
        arma::Mat<ValueType> testPoints, trialPoints;
        std::vector<ValueType> weights;

        fillDoubleSingularQuadraturePointsAndWeights(
                    desc, testPoints, trialPoints, weights);
        typedef NonseparableNumericalTestKernelTrialIntegrator<ValueType, GeometryFactory> Integrator;
        integrator = new Integrator(
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
