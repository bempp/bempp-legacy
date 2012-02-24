#include "elementary_integral_operator.hpp"

#include "assembly_options.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../common/types.hpp"
#include "../fiber/double_integrator.hpp"
#include "../fiber/integration_manager.hpp"
#include "../fiber/integration_manager_factory.hpp"
#include "../grid/entity_pointer.hpp"
#include "../grid/entity.hpp"
#include "../grid/geometry_adapter.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/grid.hpp"
#include "../space/space.hpp"

#include <armadillo>
// #include <boost/unordered_set.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <set>

namespace Bempp
{

template <typename ValueType>
std::auto_ptr<typename ElementaryIntegralOperator<ValueType>::IntegrationManager>
ElementaryIntegralOperator<ValueType>::makeIntegrationManager(
        const ElementaryIntegralOperator<ValueType>::IntegrationManagerFactory& factory) const
{
    return factory.make(testExpression(), kernel(), trialExpression());
}

template <typename ValueType>
void ElementaryIntegralOperator<ValueType>::evaluateLocalWeakForms(
        CallVariant callVariant,
        const std::vector<const EntityPointer<0>*>& elementsA,
        const EntityPointer<0>& elementB,
        LocalDofIndex localDofIndexB,
        const Space<ValueType>& spaceA,
        const Space<ValueType>& spaceB,
        ElementaryIntegralOperator<ValueType>::IntegrationManager& intMgr,
        std::vector<arma::Mat<ValueType> >& result) const
{
    typedef Fiber::DoubleIntegrator<ValueType, GeometryAdapter> Integrator;
    typedef Fiber::Basis<ValueType> Basis;

    const int elementACount = elementsA.size();
    result.resize(elementACount);

    // Get bases
    std::vector<const Basis*> basesA;
    spaceA.getBases(elementsA, basesA);
    const Basis& basisB = spaceB.basis(elementB);

    // Get index set (note: we assume that spaceA and space B refer to the same
    // grid; this has been checked in assembleWeakForm())
    std::auto_ptr<GridView> leafView(spaceA.grid().leafView());
    const IndexSet& indexSet = leafView->indexSet();

    // Construct geometry adapters (store info about geometry and vertex indices)
    boost::ptr_vector<GeometryAdapter> geometriesAOwner;
    geometriesAOwner.reserve(elementACount);
    for (int i = 0; i < elementACount; ++i)
        geometriesAOwner.push_back(
                    new GeometryAdapter(elementsA[i]->entity(), indexSet));
    std::vector<const GeometryAdapter*> geometriesA(elementACount);
    for (int i = 0; i < elementACount; ++i)
        geometriesA[i] = &geometriesAOwner[i];

    std::auto_ptr<GeometryAdapter> geometryB(
                new GeometryAdapter(elementB.entity(), indexSet));

//    std::vector<const GeometryAdapter*> geometriesA(elementACount);
//    for (int i = 0; i < elementACount; ++i)
//        geometriesA[i] = &elementsA[i]->entity().geometry();
//    const Geometry& geometryB = elementB.entity().geometry();

    // Select the integrator appropriate for each pair of elements
    std::vector<const Integrator*> integrators(elementACount);
    intMgr.getTestKernelTrialIntegrators(
                callVariant, geometriesA, *geometryB, basesA, basisB, integrators);

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and element variant (the
    // latter is important because it guarantees that all the elements have
    // the same basis functions)

    // First, find all the unique quadrature variants present
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    std::vector<QuadVariant> quadVariants(elementACount);
    for (int indexA = 0; indexA < elementACount; ++indexA)
        quadVariants[indexA] = QuadVariant(integrators[indexA], basesA[indexA]);
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<const GeometryAdapter*> activeGeometriesA;
    activeGeometriesA.reserve(elementACount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it)
    {
        const QuadVariant activeQuadVariant = *it;
        const Integrator& activeIntegrator = *it->first;
        const Basis& activeBasisA = *it->second;

        // Find all the test elements for which quadrature should proceed
        // according to the current quadrature variant
        activeGeometriesA.clear();
        for (int indexA = 0; indexA < elementACount; ++indexA)
            if (quadVariants[indexA] == activeQuadVariant)
                activeGeometriesA.push_back(geometriesA[indexA]);

        // Integrate!
        arma::Cube<ValueType> localResult;
        activeIntegrator.integrate(callVariant,
                                   activeGeometriesA, *geometryB,
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

template <typename ValueType>
void ElementaryIntegralOperator<ValueType>::evaluateLocalWeakForms(
        const std::vector<const EntityPointer<0>*>& testElements,
        const std::vector<const EntityPointer<0>*>& trialElements,
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        IntegrationManager& intMgr,
        Fiber::Array2D<arma::Mat<ValueType> >& result) const
{
    typedef Fiber::DoubleIntegrator<ValueType, GeometryAdapter> Integrator;
    typedef Fiber::Basis<ValueType> Basis;

    const int testElementCount = testElements.size();
    const int trialElementCount = trialElements.size();
    result.set_size(testElementCount, trialElementCount);

    // Get bases
    std::vector<const Basis*> testBases;
    testSpace.getBases(testElements, testBases);
    std::vector<const Basis*> trialBases;
    trialSpace.getBases(trialElements, trialBases);

    // Get index set (note: we assume that spaceA and space B refer to the same
    // grid; this has been checked in assembleWeakForm())
    std::auto_ptr<GridView> leafView(testSpace.grid().leafView());
    const IndexSet& indexSet = leafView->indexSet();

    // Construct geometry adapters (store info about geometry and vertex indices)
    boost::ptr_vector<GeometryAdapter> testGeometryOwner;
    testGeometryOwner.reserve(testElementCount);
    for (int i = 0; i < testElementCount; ++i)
        testGeometryOwner.push_back(
                    new GeometryAdapter(testElements[i]->entity(), indexSet));
    std::vector<const GeometryAdapter*> testGeometries(testElementCount);
    for (int i = 0; i < testElementCount; ++i)
        testGeometries[i] = &testGeometryOwner[i];

    // Construct geometry adapters (store info about geometry and vertex indices)
    boost::ptr_vector<GeometryAdapter> trialGeometryOwner;
    trialGeometryOwner.reserve(trialElementCount);
    for (int i = 0; i < trialElementCount; ++i)
        trialGeometryOwner.push_back(
                    new GeometryAdapter(trialElements[i]->entity(), indexSet));
    std::vector<const GeometryAdapter*> trialGeometries(trialElementCount);
    for (int i = 0; i < trialElementCount; ++i)
        trialGeometries[i] = &trialGeometryOwner[i];

    // Select the integrator appropriate for each pair of elements
    Fiber::Array2D<const Integrator*> integrators;
    intMgr.getTestKernelTrialIntegrators(
                testGeometries, trialGeometries,
                testBases, trialBases, integrators);

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and element variant (the
    // latter is important because it guarantees that all the elements have
    // the same basis functions)

    // First, find all the unique quadrature variants present
    typedef boost::tuples::tuple<const Integrator*, const Basis*, const Basis*>
            QuadVariant;
    Fiber::Array2D<QuadVariant> quadVariants(testElementCount, trialElementCount);
    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
        for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
            quadVariants(testIndex, trialIndex) =
                    QuadVariant(integrators(testIndex, trialIndex),
                                testBases[testIndex], trialBases[trialIndex]);
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    typedef std::pair<const GeometryAdapter*, const GeometryAdapter*>
            GeometryAdapterPair;
    std::vector<GeometryAdapterPair> activeGeometryPairs;
    activeGeometryPairs.reserve(testElementCount * trialElementCount);

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
        activeGeometryPairs.clear();
        for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
            for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
                if (quadVariants(testIndex, trialIndex) == activeQuadVariant)
                    activeGeometryPairs.push_back(
                                GeometryAdapterPair(testGeometries[testIndex],
                                                    trialGeometries[trialIndex]));

        // Integrate!
        arma::Cube<ValueType> localResult;
        activeIntegrator.integrate(activeGeometryPairs, activeTestBasis,
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

#ifdef COMPILE_FOR_FLOAT
template class ElementaryIntegralOperator<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class ElementaryIntegralOperator<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class ElementaryIntegralOperator<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class ElementaryIntegralOperator<std::complex<double> >;
#endif

} // namespace Bempp
