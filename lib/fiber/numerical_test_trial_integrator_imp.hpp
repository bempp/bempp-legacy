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

#include "shapeset.hpp"
#include "basis_data.hpp"
#include "conjugate.hpp"
#include "collection_of_shapeset_transformations.hpp"
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
        const CollectionOfShapesetTransformations<CoordinateType>& testTransformations,
        const CollectionOfShapesetTransformations<CoordinateType>& trialTransformations,
        const TestTrialIntegral<BasisFunctionType, ResultType>& integral,
        const OpenClHandler& openClHandler) :
    m_localQuadPoints(localQuadPoints),
    m_quadWeights(quadWeights),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testTransformations(testTransformations),
    m_trialTransformations(trialTransformations),
    m_integral(integral),
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
        const Shapeset<BasisFunctionType>& testShapeset,
        const Shapeset<BasisFunctionType>& trialShapeset,
        arma::Cube<ResultType>& result) const
{
    const size_t pointCount = m_localQuadPoints.n_cols;
    const size_t elementCount = elementIndices.size();

    // Evaluate constants
    const int componentCount = m_testTransformations.resultDimension(0);
    const int testDofCount = testShapeset.size();
    const int trialDofCount = trialShapeset.size();

    result.set_size(testDofCount, trialDofCount, elementCount);

    if (testDofCount == 0 || trialDofCount == 0 || elementCount == 0)
        return;

    if (pointCount == 0)
    {
        result.fill(0.);
        return;
    }

//    if (m_trialTransformations.codomainDimension() != componentCount)
//        throw std::runtime_error("NumericalTestTrialIntegrator::integrate(): "
//                                 "test and trial functions "
//                                 "must have the same number of components");

    BasisData<BasisFunctionType> testBasisData, trialBasisData;
    GeometricalData<CoordinateType> geomData;

    size_t testBasisDeps = 0, trialBasisDeps = 0;
    size_t geomDeps = 0; // INTEGRATION_ELEMENTS;

    m_testTransformations.addDependencies(testBasisDeps, geomDeps);
    m_trialTransformations.addDependencies(trialBasisDeps, geomDeps);
    m_integral.addGeometricalDependencies(geomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometry(m_geometryFactory.make());

    CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;

    testShapeset.evaluate(testBasisDeps, m_localQuadPoints, ALL_DOFS, testBasisData);
    trialShapeset.evaluate(trialBasisDeps, m_localQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (size_t e = 0; e < elementCount; ++e)
    {
        m_rawGeometry.setupGeometry(elementIndices[e], *geometry);
        geometry->getData(geomDeps, m_localQuadPoints, geomData);
        m_testTransformations.evaluate(testBasisData, geomData, testValues);
        m_trialTransformations.evaluate(trialBasisData, geomData, trialValues);

        m_integral.evaluate(geomData, testValues, trialValues,
                m_quadWeights, result.slice(e));

        // for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
        //     for (int testDof = 0; testDof < testDofCount; ++testDof)
        //     {
        //         ResultType sum = 0.;
        //         for (size_t point = 0; point < pointCount; ++point)
        //             for (int dim = 0; dim < componentCount; ++dim)
        //                 sum +=  m_quadWeights[point] *
        //                         geomData.integrationElements(point) *
        //                         conjugate(testValues[0](dim, testDof, point)) *
        //                         trialValues[0](dim, trialDof, point);
        //         result(testDof, trialDof, e) = sum;
        //     }
    }
}

} // namespace Fiber
