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

namespace Fiber {

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
NumericalTestTrialIntegrator<BasisFunctionType, ResultType, GeometryFactory>::
    NumericalTestTrialIntegrator(
        const Matrix<CoordinateType> &localQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const GeometryFactory &geometryFactory,
        const RawGridGeometry<CoordinateType> &rawGeometry,
        const CollectionOfShapesetTransformations<CoordinateType> &
            testTransformations,
        const CollectionOfShapesetTransformations<CoordinateType> &
            trialTransformations,
        const TestTrialIntegral<BasisFunctionType, ResultType> &integral,
        const OpenClHandler &openClHandler)
    : m_localQuadPoints(localQuadPoints), m_quadWeights(quadWeights),
      m_geometryFactory(geometryFactory), m_rawGeometry(rawGeometry),
      m_testTransformations(testTransformations),
      m_trialTransformations(trialTransformations), m_integral(integral),
      m_openClHandler(openClHandler) {
  if (localQuadPoints.cols() != quadWeights.size())
    throw std::invalid_argument("NumericalTestTrialIntegrator::"
                                "NumericalTestTrialIntegrator(): "
                                "numbers of points and weights do not match");
}

template <typename BasisFunctionType, typename ResultType,
          typename GeometryFactory>
void
NumericalTestTrialIntegrator<BasisFunctionType, ResultType, GeometryFactory>::
    integrate(const std::vector<int> &elementIndices,
              const Shapeset<BasisFunctionType> &testShapeset,
              const Shapeset<BasisFunctionType> &trialShapeset,
              std::vector<Matrix<ResultType>> &result) const {
  const size_t pointCount = m_localQuadPoints.cols();
  const size_t elementCount = elementIndices.size();

  if (pointCount == 0 || elementCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // elementCount != 0, set elements of result to 0.

  // Evaluate constants
  const int testDofCount = testShapeset.size();
  const int trialDofCount = trialShapeset.size();

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  GeometricalData<CoordinateType> geomData;

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t geomDeps = 0; // INTEGRATION_ELEMENTS;

  m_testTransformations.addDependencies(testBasisDeps, geomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, geomDeps);
  m_integral.addGeometricalDependencies(geomDeps);

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> geometry(m_geometryFactory.make());

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;

  //result.set_size(testDofCount, trialDofCount, elementCount);
  result.resize(elementCount);

  testShapeset.evaluate(testBasisDeps, m_localQuadPoints, ALL_DOFS,
                        testBasisData);
  trialShapeset.evaluate(trialBasisDeps, m_localQuadPoints, ALL_DOFS,
                         trialBasisData);

  // Iterate over the elements
  for (size_t e = 0; e < elementCount; ++e) {
    result[e].resize(testDofCount, trialDofCount);
    const int elementIndex = elementIndices[e];
    m_rawGeometry.setupGeometry(elementIndex, *geometry);
    geometry->getData(geomDeps, m_localQuadPoints, geomData);
    if (geomDeps & DOMAIN_INDEX)
      geomData.domainIndex = m_rawGeometry.domainIndex(elementIndex);
    m_testTransformations.evaluate(testBasisData, geomData, testValues);
    m_trialTransformations.evaluate(trialBasisData, geomData, trialValues);

    m_integral.evaluate(geomData, testValues, trialValues, m_quadWeights,
                        result[e]);
  }
}

} // namespace Fiber
