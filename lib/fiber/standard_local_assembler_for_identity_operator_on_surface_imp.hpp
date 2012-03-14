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
void StandardLocalAssemblerForIdentityOperatorOnSurface<ValueType, GeometryFactory>::
evaluateLocalWeakForms(
        const std::vector<int>& elementIndices,
        std::vector<arma::Mat<ValueType> >& result)
{
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
                result[e] = localResult.slice(i++);
    }
}

} // namespace Fiber
