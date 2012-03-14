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

#include "numerical_test_trial_integrator.hpp" // To keep IDEs happy

#include "basis.hpp"
#include "basis_data.hpp"
#include "expression.hpp"
#include "geometrical_data.hpp"
#include "opencl_options.hpp"
#include "raw_grid_geometry.hpp"
#include "types.hpp"

#include <cassert>
#include <memory>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
NumericalTestTrialIntegrator<ValueType, GeometryFactory>::
NumericalTestTrialIntegrator(
        const arma::Mat<ValueType>& localQuadPoints,
        const std::vector<ValueType> quadWeights,
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<ValueType>& rawGeometry,
        const Expression<ValueType>& testExpression,
        const Expression<ValueType>& trialExpression,
        const OpenClOptions& openClOptions) :
    m_localQuadPoints(localQuadPoints),
    m_quadWeights(quadWeights),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testExpression(testExpression),
    m_trialExpression(trialExpression),
    m_openClOptions(openClOptions)
{
    if (localQuadPoints.n_cols != quadWeights.size())
        throw std::invalid_argument("NumericalTestTrialIntegrator::"
                                    "NumericalTestTrialIntegrator(): "
                                    "numbers of points and weights do not match");
}

template <typename ValueType, typename GeometryFactory>
void NumericalTestTrialIntegrator<ValueType, GeometryFactory>::integrate(
        const std::vector<int>& elementIndices,
        const Basis<ValueType>& testBasis,
        const Basis<ValueType>& trialBasis,
        arma::Cube<ValueType>& result) const
{
    const int pointCount = m_localQuadPoints.n_cols;
    const int elementCount = elementIndices.size();

    if (pointCount == 0 || elementCount == 0)
        return;
    // TODO: in the (pathological) case that pointCount == 0 but
    // elementCount != 0, set elements of result to 0.

    // Evaluate constants
    const int testComponentCount = m_testExpression.codomainDimension();
    const int trialComponentCount = m_trialExpression.codomainDimension();
    const int testDofCount = testBasis.size();
    const int trialDofCount = trialBasis.size();

    BasisData<ValueType> testBasisData, trialBasisData;
    GeometricalData<ValueType> geomData;

    int testBasisDeps = 0, trialBasisDeps = 0;
    int geomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, geomDeps);
    m_trialExpression.addDependencies(trialBasisDeps, geomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometry(m_geometryFactory.make());

    arma::Cube<ValueType> testValues, trialValues;

    result.set_size(testDofCount, trialDofCount, elementCount);

    testBasis.evaluate(testBasisDeps, m_localQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (int e = 0; e < elementCount; ++e)
    {
        m_rawGeometry.setupGeometry(elementIndices[e], *geometry);
        geometry->getData(geomDeps, m_localQuadPoints, geomData);
        m_testExpression.evaluate(testBasisData, geomData, testValues);
        m_trialExpression.evaluate(trialBasisData, geomData, trialValues);

        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
            for (int testDof = 0; testDof < testDofCount; ++testDof)
            {
                ValueType sum = 0.;
                for (int point = 0; point < pointCount; ++point)
                    for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                        for (int testDim = 0; testDim < testComponentCount; ++testDim)
                            sum +=  m_quadWeights[point] *
                                    geomData.integrationElements(point) *
                                    testValues(testDim, testDof, point) *
                                    trialValues(trialDim, trialDof, point);
                result(testDof, trialDof, e) = sum;
            }
    }
}

} // namespace Fiber
