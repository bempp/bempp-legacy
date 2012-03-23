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

#include "numerical_test_function_integrator.hpp" // To keep IDEs happy

#include "basis.hpp"
#include "basis_data.hpp"
#include "expression.hpp"
#include "function.hpp"
#include "geometrical_data.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "types.hpp"

#include <stdexcept>
#include <memory>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
NumericalTestFunctionIntegrator<ValueType, GeometryFactory>::
NumericalTestFunctionIntegrator(
        const arma::Mat<ValueType>& localQuadPoints,
        const std::vector<ValueType> quadWeights,
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<ValueType>& rawGeometry,
        const Expression<ValueType>& testExpression,
        const Function<ValueType>& function,
        const OpenClHandler<ValueType,int>& openClHandler) :
    m_localQuadPoints(localQuadPoints),
    m_quadWeights(quadWeights),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testExpression(testExpression),
    m_function(function),
    m_openClHandler(openClHandler)
{
    if (localQuadPoints.n_cols != quadWeights.size())
        throw std::invalid_argument("NumericalTestTrialIntegrator::"
                                    "NumericalTestTrialIntegrator(): "
                                    "numbers of points and weights do not match");
}

template <typename ValueType, typename GeometryFactory>
void NumericalTestFunctionIntegrator<ValueType, GeometryFactory>::integrate(
        const std::vector<int>& elementIndices,
        const Basis<ValueType>& testBasis,
        arma::Mat<ValueType>& result) const
{
    const int pointCount = m_localQuadPoints.n_cols;
    const int elementCount = elementIndices.size();

    if (pointCount == 0 || elementCount == 0)
        return;
    // TODO: in the (pathological) case that pointCount == 0 but
    // elementCount != 0, set elements of result to 0.

    // Evaluate constants
    const int componentCount = m_testExpression.codomainDimension();
    const int testDofCount = testBasis.size();

    if (m_function.codomainDimension() != componentCount)
        throw std::runtime_error("NumericalTestFunctionIntegrator::integrate(): "
                                 "test functions and the \"arbitrary\" function "
                                 "must have the same number of components");

    BasisData<ValueType> testBasisData;
    GeometricalData<ValueType> geomData;

    int testBasisDeps = 0;
    int geomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, geomDeps);
    m_function.addGeometricalDependencies(geomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometry(m_geometryFactory.make());

    arma::Cube<ValueType> testValues;
    arma::Mat<ValueType> functionValues;

    result.set_size(testDofCount, elementCount);

    testBasis.evaluate(testBasisDeps, m_localQuadPoints, ALL_DOFS, testBasisData);

    // Iterate over the elements
    for (int e = 0; e < elementCount; ++e)
    {
        m_rawGeometry.setupGeometry(elementIndices[e], *geometry);
        geometry->getData(geomDeps, m_localQuadPoints, geomData);
        m_testExpression.evaluate(testBasisData, geomData, testValues);
        m_function.evaluate(geomData, functionValues);

        for (int testDof = 0; testDof < testDofCount; ++testDof)
        {
            ValueType sum = 0.;
            for (int point = 0; point < pointCount; ++point)
                for (int dim = 0; dim < componentCount; ++dim)
                    sum +=  m_quadWeights[point] *
                            geomData.integrationElements(point) *
                            testValues(dim, testDof, point) *
                            functionValues(dim, point);
            result(testDof, e) = sum;
        }
    }
}

} // namespace Fiber
