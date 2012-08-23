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

// Keep IDEs happy
#include "default_local_assembler_for_grid_functions_on_surfaces.hpp"

#include "../common/common.hpp"

#include "numerical_test_function_integrator.hpp"

#include <set>
#include <utility>

namespace Fiber
{

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
DefaultLocalAssemblerForGridFunctionsOnSurfaces<
BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>::
DefaultLocalAssemblerForGridFunctionsOnSurfaces(
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& testBases,
        const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& testTransformations,
        const shared_ptr<const Function<UserFunctionType> >& function,
        const shared_ptr<const OpenClHandler>& openClHandler) :
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testBases(testBases),
    m_testTransformations(testTransformations),
    m_function(function),
    m_openClHandler(openClHandler)
{
    if (rawGeometry->vertices().n_rows != 3)
        throw std::invalid_argument(
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces::"
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces(): "
                "vertex coordinates must be three-dimensional");
    if (rawGeometry->elementCornerIndices().n_rows < 3 ||
            4 < rawGeometry->elementCornerIndices().n_rows)
        throw std::invalid_argument(
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces::"
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces(): "
                "all elements must be triangular or quadrilateral");
    const size_t elementCount = rawGeometry->elementCornerIndices().n_cols;
    if (!rawGeometry->auxData().is_empty() &&
            rawGeometry->auxData().n_cols != elementCount)
        throw std::invalid_argument(
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces::"
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces(): "
                "number of columns of auxData must match that of "
                "elementCornerIndices");
    if (testBases->size() != elementCount)
        throw std::invalid_argument(
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces::"
                "DefaultLocalAssemblerForGridFunctionsOnSurfaces(): "
                "size of testBases must match the number of columns of "
                "elementCornerIndices");
}

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
DefaultLocalAssemblerForGridFunctionsOnSurfaces<
BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>::
~DefaultLocalAssemblerForGridFunctionsOnSurfaces()
{
    // Note: obviously the destructor is assumed to be called only after
    // all threads have ceased using the assembler!

    for (typename IntegratorMap::const_iterator it = m_testFunctionIntegrators.begin();
         it != m_testFunctionIntegrators.end(); ++it)
        delete it->second;
    m_testFunctionIntegrators.clear();
}

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForGridFunctionsOnSurfaces<
BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& elementIndices,
        std::vector<arma::Col<ResultType> >& result)
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
                QuadVariant(integrator, (*m_testBases)[activeTestElementIndex]);
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
    }
}

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
const TestFunctionIntegrator<BasisFunctionType, ResultType>&
DefaultLocalAssemblerForGridFunctionsOnSurfaces<
BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>::
selectIntegrator(int elementIndex)
{
    SingleQuadratureDescriptor desc;

    // Get number of corners of the specified element
    desc.vertexCount = m_rawGeometry->elementCornerCount(elementIndex);

    // Determine integrand's order and required quadrature order
    desc.order = (*m_testBases)[elementIndex]->order() +
            orderIncrement(elementIndex);

    return getIntegrator(desc);
}

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
const TestFunctionIntegrator<BasisFunctionType, ResultType>&
DefaultLocalAssemblerForGridFunctionsOnSurfaces<
BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>::
getIntegrator(const SingleQuadratureDescriptor& desc)
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
    fillSingleQuadraturePointsAndWeights(desc.vertexCount, desc.order,
                                         points, weights);

    typedef NumericalTestFunctionIntegrator<BasisFunctionType, UserFunctionType,
            ResultType, GeometryFactory> ConcreteIntegrator;
    Integrator* integrator(
                new ConcreteIntegrator(points, weights,
                                       *m_geometryFactory, *m_rawGeometry,
                                       *m_testTransformations, *m_function,
                                       *m_openClHandler));

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

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
inline int
DefaultLocalAssemblerForGridFunctionsOnSurfaces<
BasisFunctionType, UserFunctionType, ResultType, GeometryFactory>::
orderIncrement(int elementIndex) const
{
    // TODO: add to constructor an option for increased-order quadrature
    return (*m_testBases)[elementIndex]->order() + 1;
}

} // namespace Fiber
