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

#include "numerical_test_function_integrator.hpp" // To keep IDEs happy

#include "../common/common.hpp"

#include "shapeset.hpp"
#include "basis_data.hpp"
#include "collection_of_shapeset_transformations.hpp"
#include "conjugate.hpp"
#include "function.hpp"
#include "geometrical_data.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "types.hpp"

#include <stdexcept>
#include <memory>

namespace Fiber {

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
NumericalTestFunctionIntegrator<BasisFunctionType, UserFunctionType, ResultType,
                                GeometryFactory>::
    NumericalTestFunctionIntegrator(
        const Matrix<CoordinateType> &localQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const GeometryFactory &geometryFactory,
        const RawGridGeometry<CoordinateType> &rawGeometry,
        const CollectionOfShapesetTransformations<CoordinateType> &
            testTransformations,
        const Function<UserFunctionType> &function,
        const OpenClHandler &openClHandler)
    : m_localQuadPoints(localQuadPoints), m_quadWeights(quadWeights),
      m_geometryFactory(geometryFactory), m_rawGeometry(rawGeometry),
      m_testTransformations(testTransformations), m_function(function),
      m_openClHandler(openClHandler) {
  if (localQuadPoints.cols() != quadWeights.size())
    throw std::invalid_argument("NumericalTestTrialIntegrator::"
                                "NumericalTestTrialIntegrator(): "
                                "numbers of points and weights do not match");

  if (testTransformations.transformationCount() != 1)
    throw std::invalid_argument("NumericalTestTrialIntegrator::"
                                "NumericalTestTrialIntegrator(): "
                                "test transformation collection "
                                "must contain exactly one element");
}

template <typename BasisFunctionType, typename UserFunctionType,
          typename ResultType, typename GeometryFactory>
void NumericalTestFunctionIntegrator<
    BasisFunctionType, UserFunctionType, ResultType,
    GeometryFactory>::integrate(const std::vector<int> &elementIndices,
                                const Shapeset<BasisFunctionType> &testShapeset,
                                Matrix<ResultType> &result) const {
  const size_t pointCount = m_localQuadPoints.cols();
  const size_t elementCount = elementIndices.size();

  if (pointCount == 0 || elementCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // elementCount != 0, set elements of result to 0.

  // Evaluate constants
  const int componentCount = m_testTransformations.resultDimension(0);
  const int testDofCount = testShapeset.size();

  if (m_function.codomainDimension() != componentCount)
    throw std::runtime_error("NumericalTestFunctionIntegrator::integrate(): "
                             "test functions and the \"arbitrary\" function "
                             "must have the same number of components");

  BasisData<BasisFunctionType> testBasisData;
  GeometricalData<CoordinateType> geomData;

  size_t testBasisDeps = 0;
  size_t geomDeps = INTEGRATION_ELEMENTS;

  m_testTransformations.addDependencies(testBasisDeps, geomDeps);
  m_function.addGeometricalDependencies(geomDeps);

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> geometry(m_geometryFactory.make());

  Fiber::CollectionOf3dArrays<BasisFunctionType> testValues;
  Matrix<UserFunctionType> functionValues;

  result.resize(testDofCount, elementCount);

  testShapeset.evaluate(testBasisDeps, m_localQuadPoints, ALL_DOFS,
                        testBasisData);

  // Iterate over the elements
  for (size_t e = 0; e < elementCount; ++e) {
    const int elementIndex = elementIndices[e];
    m_rawGeometry.setupGeometry(elementIndex, *geometry);
    geometry->getData(geomDeps, m_localQuadPoints, geomData);
    if (geomDeps & DOMAIN_INDEX)
      geomData.domainIndex = m_rawGeometry.domainIndex(elementIndex);
    m_testTransformations.evaluate(testBasisData, geomData, testValues);
    m_function.evaluate(geomData, functionValues);

    for (int testDof = 0; testDof < testDofCount; ++testDof) {
      ResultType sum = 0.;
      for (size_t point = 0; point < pointCount; ++point)
        for (int dim = 0; dim < componentCount; ++dim)
          sum += m_quadWeights[point] * geomData.integrationElements(point) *
                 conjugate(testValues[0](dim, testDof, point)) *
                 functionValues(dim, point);
      result(testDof, e) = sum;
    }
  }
}

} // namespace Fiber
