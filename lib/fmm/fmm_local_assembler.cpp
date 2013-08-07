
#include "fmm_local_assembler.hpp"

#include "../space/space.hpp"
#include "../assembly/assembly_options.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../assembly/local_assembler_construction_helper.hpp"
#include "../fiber/function.hpp"

#include "../fiber/collection_of_3d_arrays.hpp"
#include "../fiber/numerical_test_function_integrator.hpp"
#include "../fiber/numerical_trial_function_integrator.hpp"
#include "../fiber/test_function_integrator.hpp"
#include "../fiber/quadrature_options.hpp"
#include "../grid/geometry_factory.hpp" // defines GeometryFactory

#include <boost/type_traits/is_complex.hpp>

#include <set>

namespace Bempp
{

/*
 * Constructor inspired by the following, since they geometry etc are not
 * obtained in the usual way from the elementary operator
 * /assembly/grid_function.cpp
 * /fiber/default_local_assembler_for_grid_functions_on_surfaces_imp.hpp
 * /assembly/abstract_boundary_operator.cpp
 * constuct once per process per test or trial space
 */
template <typename BasisFunctionType, typename ResultType>
FmmLocalAssembler<BasisFunctionType, ResultType>::FmmLocalAssembler(
	const Space<BasisFunctionType>& space, 
	const AssemblyOptions& options, bool conjugateBasis)
	: m_conjugateBasis(conjugateBasis)
{
	// ADD some initialisation checks here

	// initialise objects required to make integrator
	// from /assembly/grid_function.cpp
	typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
	typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;
	typedef LocalAssemblerConstructionHelper Helper;

	// Collect grid data, formerly this worked on the dualToRange space = test
	Helper::collectGridData(*space.grid(),
		m_rawGeometry, m_geometryFactory);

	// Get reference to the basis transformation
	m_transformations = make_shared_from_ref(space.shapeFunctionValue());

	Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                              m_rawGeometry, m_openClHandler);

	Helper::collectBases(space, m_bases);
}

/*
 * The remainder follows the form of 
 * default_local_assembler_for_grid_functions_on_surfaces_imp.hpp
 * with the expection that we have loops over the spherical quadrature points
 */

template <typename BasisFunctionType, typename ResultType>
FmmLocalAssembler<BasisFunctionType, ResultType>::~FmmLocalAssembler()
{
	clearIntegratorMap();
}

template <typename BasisFunctionType, typename ResultType>
void FmmLocalAssembler<BasisFunctionType, ResultType>::clearIntegratorMap()
{
    // Note: obviously the destructor is assumed to be called only after
    // all threads have ceased using the assembler!

	if (!m_testFunctionIntegrators.empty())
	{
		for (typename IntegratorMap::const_iterator 
			it = m_testFunctionIntegrators.begin();
			it != m_testFunctionIntegrators.end(); ++it)
			delete it->second;
		m_testFunctionIntegrators.clear();
	}
}


template <typename BasisFunctionType, typename ResultType>
const Fiber::TestFunctionIntegrator<BasisFunctionType, ResultType>& 
FmmLocalAssembler<BasisFunctionType, ResultType>::selectIntegrator(int elementIndex)
{
	Fiber::SingleQuadratureDescriptor desc;

	// Get number of corners of the specified element
	desc.vertexCount = m_rawGeometry->elementCornerCount(elementIndex);

	const int defaultOrder = 2 * (*m_bases)[elementIndex]->order() + 1;
	Fiber::QuadratureOptions quadratureOptions;		// use default quadrature (where is user var?)
	quadratureOptions.setRelativeQuadratureOrder(2);	// uncomment to increase accuracy
	desc.order = quadratureOptions.quadratureOrder(defaultOrder);

	return getIntegrator(desc);
}

template <typename BasisFunctionType, typename ResultType>
const Fiber::TestFunctionIntegrator<BasisFunctionType, ResultType>& 
FmmLocalAssembler<BasisFunctionType, ResultType>::getIntegrator(
	const Fiber::SingleQuadratureDescriptor& desc)
{
	typename IntegratorMap::iterator it = m_testFunctionIntegrators.find(desc);
	if (it != m_testFunctionIntegrators.end())
	{
		// std::cout << "getIntegrator(: " << index << "): integrator found" << std::endl;
		return *it->second;
	}
	// std::cout << "getIntegrator(: " << index << "): integrator not found" << std::endl;

	// Integrator doesn't exist yet and must be created.
	arma::Mat<CoordinateType> points;
	std::vector<CoordinateType> weights;

	Fiber::fillSingleQuadraturePointsAndWeights(desc.vertexCount, desc.order,
		                               points, weights);


	// don't currently know the concrete type for Geometry factory, can improve speed if known
	Integrator* integrator;

	if (m_conjugateBasis) {
		typedef Fiber::NumericalTestFunctionIntegrator<BasisFunctionType, UserFunctionType,
			ResultType, GeometryFactory> ConcreteIntegrator;

			integrator = new ConcreteIntegrator(points, weights,
                                       *m_geometryFactory, *m_rawGeometry,
                                       *m_transformations, *m_function,
                                       *m_openClHandler );
	} else {
		typedef Fiber::NumericalTrialFunctionIntegrator<BasisFunctionType, UserFunctionType,
			ResultType, GeometryFactory> ConcreteIntegrator;

			integrator = new ConcreteIntegrator(points, weights,
                                       *m_geometryFactory, *m_rawGeometry,
                                       *m_transformations, *m_function,
                                       *m_openClHandler );
	}

	// Attempt to insert the newly created integrator into the map
	std::pair<typename IntegratorMap::iterator, bool> result =
		m_testFunctionIntegrators.insert(std::make_pair(desc, integrator));
	if (result.second)
		// Insertion succeeded. The newly created integrator will be deleted in
		// our own destructor
		;
	else
		// Insertion failed -- another thread was faster. Delete the newly
		// created integrator.
		delete integrator;

	// Return pointer to the integrator that ended up in the map.
	return *result.first->second;
}

// modified to return a vector of results, one for each quadrature point on sphere
template <typename BasisFunctionType, typename ResultType>
void FmmLocalAssembler<BasisFunctionType, ResultType>::evaluateLocalWeakForms(
	const std::vector<int>& elementIndices,
	std::vector<arma::Col<ResultType> > & result)
{
	typedef Fiber::Basis<BasisFunctionType> Basis;

	const int elementCount = elementIndices.size();
	result.resize(elementCount);

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    std::vector<QuadVariant> quadVariants(elementCount);

	for (int testIndex = 0; testIndex < elementCount; ++testIndex)
	{
		const int activeTestElementIndex = elementIndices[testIndex];
		const Integrator* integrator =
			&selectIntegrator(activeTestElementIndex);
		quadVariants[testIndex] =
			QuadVariant(integrator, (*m_bases)[activeTestElementIndex]);
	}

	// Integration will proceed in batches of element pairs having the same
	// "quadrature variant", i.e. integrator and test basis

	// Find all the unique quadrature variants present
	typedef std::set<QuadVariant> QuadVariantSet;
	// Set of unique quadrature variants
	QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

	std::vector<int> activeElementIndices;
	activeElementIndices.reserve(elementCount);

    // Now loop over unique quadrature variants
	for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
		it != uniqueQuadVariants.end(); ++it)
	{
		const QuadVariant activeQuadVariant = *it;
		const Integrator& activeIntegrator = *it->first;
		const Basis& activeTestBasis = *it->second;

		// Find all the test elements for which quadrature should proceed
		// according to the current quadrature variant
		activeElementIndices.clear();
		for (int e = 0; e < elementCount; ++e)
		if (quadVariants[e] == activeQuadVariant)
			activeElementIndices.push_back(elementIndices[e]);

		// Integrate!
		arma::Mat<ResultType> localResult;
		activeIntegrator.integrate(activeElementIndices,
					activeTestBasis,
					localResult);

		// Distribute the just calculated integrals into the result array
		// that will be returned to caller
		int i = 0;
		for (int e = 0; e < elementCount; ++e)
			if (quadVariants[e] == activeQuadVariant)
				result[e] = localResult.col(i++);
	} // for each quadrature variant
}


template <typename BasisFunctionType, typename ResultType>
void FmmLocalAssembler<BasisFunctionType, ResultType>::setFunction(
	Fiber::Function<ResultType> *function)
{
	clearIntegratorMap();
	m_function = function;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmLocalAssembler);

} // namespace Bempp
