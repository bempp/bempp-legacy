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
#include "../space/space.hpp"

#include <armadillo>
// #include <boost/unordered_set.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <set>

namespace Bempp
{

template <typename ValueType>
std::auto_ptr<Fiber::IntegrationManager<ValueType, Geometry> >
ElementaryIntegralOperator<ValueType>::makeIntegrationManager(
        const Fiber::IntegrationManagerFactory<ValueType, Geometry>& factory) const
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
        Fiber::IntegrationManager<ValueType, Geometry>& integrationMgr,
        std::vector<arma::Mat<ValueType> >& result) const
{
    typedef Fiber::DoubleIntegrator<ValueType, Geometry> Integrator;
    typedef Fiber::Basis<ValueType> Basis;

    const int elementACount = elementsA.size();
    result.resize(elementACount);

    // Get bases
    std::vector<const Basis*> basesA;
    spaceA.getBases(elementsA, basesA);
    const Basis& basisB = spaceB.basis(elementB);

    // Get geometries
    std::vector<const Geometry*> geometriesA(elementACount);
    for (int i = 0; i < elementACount; ++i)
        geometriesA[i] = &elementsA[i]->entity().geometry();
    const Geometry& geometryB = elementB.entity().geometry();

    // Select the integrator appropriate for each pair of elements
    std::vector<const Integrator*> integrators(elementACount);
    integrationMgr.getTestKernelTrialIntegrators(
                callVariant, geometriesA, geometryB, basesA, basisB, integrators);

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and element variant (the
    // latter is important because it guarantees that all the elements have
    // the same basis functions)

    // First, find all the unique quadrature variants present
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    std::vector<QuadVariant> quadVariants(elementACount); // Temporary vector
    for (int indexA = 0; indexA < elementACount; ++indexA)
        quadVariants[indexA] = QuadVariant(integrators[indexA], basesA[indexA]);
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<const Geometry*> activeGeometriesA;
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
                                   activeGeometriesA, geometryB,
                                   activeBasisA, basisB, ALL_DOFS,
                                   localResult);

        // Distribute the just calculated integrals into the result array
        // that will be returned to caller
        int i = 0;
        for (int indexA = 0; indexA < elementACount; ++indexA)
            if (quadVariants[indexA] == activeQuadVariant)
                result[indexA] = localResult.slice(i++);
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
