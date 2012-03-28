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
#include "standard_local_assembler_for_identity_operator_on_surface.hpp"

#include <boost/tuple/tuple_comparison.hpp>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
StandardLocalAssemblerForIdentityOperatorOnSurface(
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<ValueType>& rawGeometry,
        const std::vector<const Basis<ValueType>*>& testBases,
        const std::vector<const Basis<ValueType>*>& trialBases,
        const Expression<ValueType>& testExpression,
        const Expression<ValueType>& trialExpression,
        ValueType multiplier,
        const OpenClHandler<ValueType,int>& openClHandler) :
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testBases(testBases),
    m_trialBases(trialBases),
    m_testExpression(testExpression),
    m_trialExpression(trialExpression),
    m_multiplier(multiplier),
    m_openClHandler(openClHandler)
{
    if (rawGeometry.vertices().n_rows != 3)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                "vertex coordinates must be three-dimensional");
    const int elementCount = rawGeometry.elementCornerIndices().n_cols;
    if (rawGeometry.elementCornerIndices().n_rows < 3 ||
            4 < rawGeometry.elementCornerIndices().n_rows)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                "Elements must have either 3 or 4 corners");
    if (!rawGeometry.auxData().is_empty() &&
            rawGeometry.auxData().n_cols != elementCount)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                "number of columns of auxData must match that of "
                "elementCornerIndices");
    if (testBases.size() != elementCount)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                "size of testBases must match the number of columns of "
                "elementCornerIndices");
    if (trialBases.size() != elementCount)
        throw std::invalid_argument(
                "StandardLocalAssemblerForIdentityOperatorOnSurface::"
                "StandardLocalAssemblerForIdentityOperatorOnSurface(): "
                "size of trialBases must match the number of columns of "
                "elementCornerIndices");
}

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
evaluateLocalWeakForms(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        LocalDofIndex localDofIndexB,
        std::vector<arma::Mat<ValueType> >& result)
{
    // Probably will never be called
    throw std::runtime_error("StandardLocalAssemblerForIdentityOperatorOnSurface::"
                             "evaluateLocalWeakForms(): "
                             "this overload not implemented yet");
}

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& testElementIndices,
        const std::vector<int>& trialElementIndices,
        Fiber::Array2D<arma::Mat<ValueType> >& result)
{
    // Probably will never be called
    throw std::runtime_error("StandardLocalAssemblerForIdentityOperatorOnSurface::"
                             "evaluateLocalWeakForms(): "
                             "this overload not implemented yet");
}

template <typename ValueType, typename GeometryFactory>
void
StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& elementIndices,
        std::vector<arma::Mat<ValueType> >& result)
{
    // The only overload likely to be needed for identity operators
    typedef TestTrialIntegrator<ValueType> Integrator;
    typedef Basis<ValueType> Basis;

    const int elementCount = elementIndices.size();
    result.resize(elementCount);

    // Select integrator for each element
    typedef boost::tuples::tuple<const Integrator*, const Basis*, const Basis*>
            QuadVariant;
    std::vector<QuadVariant> quadVariants(elementCount);
    for (int i = 0; i < elementCount; ++i)
    {
        const Integrator* integrator = &selectIntegrator(elementIndices[i]);
        quadVariants[i] = QuadVariant(integrator, m_testBases[elementIndices[i]],
                                      m_trialBases[elementIndices[i]]);
    }

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and bases

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
        const Integrator& activeIntegrator = *it->template get<0>();
        const Basis& activeTestBasis  = *it->template get<1>();
        const Basis& activeTrialBasis = *it->template get<2>();

        // Find all the test elements for which quadrature should proceed
        // according to the current quadrature variant
        activeElementIndices.clear();
        for (int e = 0; e < elementCount; ++e)
            if (quadVariants[e] == activeQuadVariant)
                activeElementIndices.push_back(elementIndices[e]);

        // Integrate!
        arma::Cube<ValueType> localResult;
        activeIntegrator.integrate(activeElementIndices,
                                   activeTestBasis, activeTrialBasis,
                                   localResult);

        // Distribute the just calculated integrals into the result array
        // that will be returned to caller
        int i = 0;
        for (int e = 0; e < elementCount; ++e)
            if (quadVariants[e] == activeQuadVariant)
                result[e] = m_multiplier * localResult.slice(i++);
    }
}

template <typename ValueType, typename GeometryFactory>
const TestTrialIntegrator<ValueType>&
StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
selectIntegrator(int elementIndex)
{
    SingleQuadratureDescriptor desc;

    // Get number of corners of the specified element
    desc.vertexCount = m_rawGeometry.elementCornerCount(elementIndex);

    // Determine integrand's order and required quadrature order
    const int expressionOrder =
            m_testBases[elementIndex]->order() +
            m_trialBases[elementIndex]->order();
    desc.order = ((expressionOrder + 1) + 1 /* round up */) / 2;

    return getIntegrator(desc);
}

template <typename ValueType, typename GeometryFactory>
const TestTrialIntegrator<ValueType>&
StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
getIntegrator(const SingleQuadratureDescriptor& desc)
{
    typename IntegratorMap::iterator it = m_testTrialIntegrators.find(desc);
    if (it != m_testTrialIntegrators.end())
    {
//            std::cout << "getIntegrator(: " << index << "): integrator found" << std::endl;
        return *it->second;
    }
//        std::cout << "getIntegrator(: " << index << "): integrator not found" << std::endl;

    // Integrator doesn't exist yet and must be created.
    arma::Mat<ValueType> points;
    std::vector<ValueType> weights;
    fillSingleQuadraturePointsAndWeights(desc.vertexCount, desc.order,
                                         points, weights);

    typedef NumericalTestTrialIntegrator<ValueType, GeometryFactory> Integrator;
    std::auto_ptr<TestTrialIntegrator<ValueType> > integrator(
            new Integrator(points, weights,
                           m_geometryFactory, m_rawGeometry,
                           m_testExpression, m_trialExpression,
                           m_openClHandler));

    return *m_testTrialIntegrators.insert(desc, integrator).first->second;
}

} // namespace Fiber
