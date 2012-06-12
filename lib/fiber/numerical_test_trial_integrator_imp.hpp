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

#include "numerical_test_trial_integrator.hpp" // To keep IDEs happy

#include "basis.hpp"
#include "basis_data.hpp"
#include "conjugate.hpp"
#include "expression.hpp"
#include "geometrical_data.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "types.hpp"

#include <cassert>
#include <memory>

namespace Fiber
{

template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
NumericalTestTrialIntegrator<BasisFunctionType, ResultType, GeometryFactory>::
NumericalTestTrialIntegrator(
        const arma::Mat<CoordinateType>& localQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const Expression<CoordinateType>& testExpression,
        const Expression<CoordinateType>& trialExpression,
        const OpenClHandler& openClHandler) :
    m_localQuadPoints(localQuadPoints),
    m_quadWeights(quadWeights),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testExpression(testExpression),
    m_trialExpression(trialExpression),
    m_openClHandler(openClHandler)
{
    if (localQuadPoints.n_cols != quadWeights.size())
        throw std::invalid_argument("NumericalTestTrialIntegrator::"
                                    "NumericalTestTrialIntegrator(): "
                                    "numbers of points and weights do not match");
}

template <typename BasisFunctionType, typename ResultType, typename GeometryFactory>
void NumericalTestTrialIntegrator<BasisFunctionType, ResultType, GeometryFactory>::integrate(
        const std::vector<int>& elementIndices,
        const Basis<BasisFunctionType>& testBasis,
        const Basis<BasisFunctionType>& trialBasis,
        arma::Cube<ResultType>& result) const
{
    const size_t pointCount = m_localQuadPoints.n_cols;
    const size_t elementCount = elementIndices.size();

    if (pointCount == 0 || elementCount == 0)
        return;
    // TODO: in the (pathological) case that pointCount == 0 but
    // elementCount != 0, set elements of result to 0.

    // Evaluate constants
    const int componentCount = m_testExpression.codomainDimension();
    const int testDofCount = testBasis.size();
    const int trialDofCount = trialBasis.size();

    if (m_trialExpression.codomainDimension() != componentCount)
        throw std::runtime_error("NumericalTestTrialIntegrator::integrate(): "
                                 "test and trial functions "
                                 "must have the same number of components");

    BasisData<BasisFunctionType> testBasisData, trialBasisData;
    GeometricalData<CoordinateType> geomData;

    size_t testBasisDeps = 0, trialBasisDeps = 0;
    size_t geomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, geomDeps);
    m_trialExpression.addDependencies(trialBasisDeps, geomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometry(m_geometryFactory.make());

    arma::Cube<BasisFunctionType> testValues, trialValues;

    result.set_size(testDofCount, trialDofCount, elementCount);

    testBasis.evaluate(testBasisDeps, m_localQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (size_t e = 0; e < elementCount; ++e)
    {
        m_rawGeometry.setupGeometry(elementIndices[e], *geometry);
        geometry->getData(geomDeps, m_localQuadPoints, geomData);
        m_testExpression.evaluate(testBasisData, geomData, testValues);
        m_trialExpression.evaluate(trialBasisData, geomData, trialValues);

        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
            for (int testDof = 0; testDof < testDofCount; ++testDof)
            {
                ResultType sum = 0.;
                for (size_t point = 0; point < pointCount; ++point)
                    for (int dim = 0; dim < componentCount; ++dim)
                        sum +=  m_quadWeights[point] *
                                geomData.integrationElements(point) *
                                conjugate(testValues(dim, testDof, point)) *
                                trialValues(dim, trialDof, point);
                result(testDof, trialDof, e) = sum;
            }
    }
}

} // namespace Fiber
