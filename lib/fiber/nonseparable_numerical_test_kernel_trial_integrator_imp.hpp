// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "nonseparable_numerical_test_kernel_trial_integrator.hpp" // To keep IDEs happy

#include "_3d_array.hpp"
#include "_4d_array.hpp"
#include "collection_of_3d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"

#include "basis.hpp"
#include "basis_data.hpp"
#include "conjugate.hpp"
#include "collection_of_basis_transformations.hpp"
#include "geometrical_data.hpp"
#include "collection_of_kernels.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "test_kernel_trial_integral.hpp"
#include "types.hpp"

#include <cassert>
#include <iostream>
#include <memory>

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
NonseparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                               ResultType, GeometryFactory>::
    NonseparableNumericalTestKernelTrialIntegrator(
        const arma::Mat<CoordinateType> &localTestQuadPoints,
        const arma::Mat<CoordinateType> &localTrialQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const GeometryFactory &testGeometryFactory,
        const GeometryFactory &trialGeometryFactory,
        const RawGridGeometry<CoordinateType> &testRawGeometry,
        const RawGridGeometry<CoordinateType> &trialRawGeometry,
        const CollectionOfShapesetTransformations<CoordinateType> &
            testTransformations,
        const CollectionOfKernels<KernelType> &kernels,
        const CollectionOfShapesetTransformations<CoordinateType> &
            trialTransformations,
        const TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                      ResultType> &integral,
        const OpenClHandler &openClHandler)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints), m_quadWeights(quadWeights),
      m_testGeometryFactory(testGeometryFactory),
      m_trialGeometryFactory(trialGeometryFactory),
      m_testRawGeometry(testRawGeometry), m_trialRawGeometry(trialRawGeometry),
      m_testTransformations(testTransformations), m_kernels(kernels),
      m_trialTransformations(trialTransformations), m_integral(integral),
      m_openClHandler(openClHandler) {
  const size_t pointCount = quadWeights.size();
  if (localTestQuadPoints.n_cols != pointCount ||
      localTrialQuadPoints.n_cols != pointCount)
    throw std::invalid_argument(
        "NonseparableNumericalTestKernelTrialIntegrator::"
        "NonseparableNumericalTestKernelTrialIntegrator(): "
        "numbers of points and weights do not match");
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
NonseparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::~NonseparableNumericalTestKernelTrialIntegrator() {
  // Note: obviously the destructor is assumed to be called only after
  // all threads have ceased using the integrator!

  for (typename BasisDataCache::const_iterator it =
           m_cachedTestBasisData.begin();
       it != m_cachedTestBasisData.end(); ++it)
    delete it->second;
  m_cachedTestBasisData.clear();
  for (typename BasisDataCache::const_iterator it =
           m_cachedTrialBasisData.begin();
       it != m_cachedTrialBasisData.end(); ++it)
    delete it->second;
  m_cachedTrialBasisData.clear();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
const BasisData<BasisFunctionType> &
NonseparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::basisData(ElementType type,
                                const Shapeset<BasisFunctionType> &shapeset)
    const {
  BasisDataCache &cache =
      type == TEST ? m_cachedTestBasisData : m_cachedTrialBasisData;
  const CollectionOfShapesetTransformations<CoordinateType> &transformations =
      type == TEST ? m_testTransformations : m_trialTransformations;
  const arma::Mat<CoordinateType> &localQuadPoints =
      type == TEST ? m_localTestQuadPoints : m_localTrialQuadPoints;

  typename BasisDataCache::iterator it = cache.find(&shapeset);
  if (it == cache.end()) {
    // FIXME: At the beginning of each loop, all threads want to calculate
    // these basis data. It would be good to do singular integral
    // calculation in such a way that all types would be merged in a single
    // parallel loop; this should solve this problem.

    // Need to calculate basis data from scratch
    std::unique_ptr<BasisData<BasisFunctionType>> basisData(
        new BasisData<BasisFunctionType>);
    size_t basisDeps = 0, geomDeps = 0;
    transformations.addDependencies(basisDeps, geomDeps);
    shapeset.evaluate(basisDeps, localQuadPoints, ALL_DOFS, *basisData);
    // Release the newly created basis data object from the auto pointer and
    // // try to insert it into the map. If the insertion succeeds, the map
    // // takes ownership of the object, which will be deleted in the
    // // destructor of the integrator. If the insertion fails, it means
    // // another thread has already inserted the basis data for this shapeset,
    // // so our basis data object is now redundant and should be released.
    BasisData<BasisFunctionType> *ptrBasisData = basisData.release();
    std::pair<typename BasisDataCache::iterator, bool> result =
        cache.insert(std::make_pair(&shapeset, ptrBasisData));
    it = result.first;  // this is the object that ended up into the cache
    if (!result.second) // insertion failed
      delete ptrBasisData;
  }
  return *it->second;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void NonseparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    integrate(CallVariant callVariant, const std::vector<int> &elementIndicesA,
              int elementIndexB, const Shapeset<BasisFunctionType> &basisA,
              const Shapeset<BasisFunctionType> &basisB,
              LocalDofIndex localDofIndexB,
              const std::vector<arma::Mat<ResultType> *> &result) const {
  const int pointCount = m_quadWeights.size();
  const int elementACount = elementIndicesA.size();

  if (result.size() != elementIndicesA.size())
    throw std::invalid_argument(
        "NonseparableNumericalTestKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'elementIndicesA' must have the same number "
        "of elements");
  if (pointCount == 0 || elementACount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryCount != 0, set elements of result to 0.

  // Evaluate constants

  const int dofCountA = basisA.size();
  const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
  const int testDofCount = callVariant == TEST_TRIAL ? dofCountA : dofCountB;
  const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB : dofCountA;

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  GeometricalData<CoordinateType> &testGeomData = m_testGeomData.local();
  GeometricalData<CoordinateType> &trialGeomData = m_trialGeomData.local();

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  typedef typename GeometryFactory::Geometry Geometry;

  std::unique_ptr<Geometry> geometryA, geometryB;
  const RawGridGeometry<CoordinateType> *rawGeometryA = 0, *rawGeometryB = 0;
  if (callVariant == TEST_TRIAL) {
    geometryA = m_testGeometryFactory.make();
    geometryB = m_trialGeometryFactory.make();
    rawGeometryA = &m_testRawGeometry;
    rawGeometryB = &m_trialRawGeometry;
  } else {
    geometryA = m_trialGeometryFactory.make();
    geometryB = m_testGeometryFactory.make();
    rawGeometryA = &m_trialRawGeometry;
    rawGeometryB = &m_testRawGeometry;
  }

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
  CollectionOf3dArrays<KernelType> kernelValues;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->set_size(testDofCount, trialDofCount);
  }

  rawGeometryB->setupGeometry(elementIndexB, *geometryB);
  if (callVariant == TEST_TRIAL) {
    basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                    testBasisData);
    basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB,
                    trialBasisData);
    geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
    if (trialGeomDeps & DOMAIN_INDEX)
      trialGeomData.domainIndex = rawGeometryB->domainIndex(elementIndexB);
    m_trialTransformations.evaluate(trialBasisData, trialGeomData, trialValues);
  } else {
    basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                    trialBasisData);
    basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB,
                    testBasisData);
    geometryB->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
    if (testGeomDeps & DOMAIN_INDEX)
      testGeomData.domainIndex = rawGeometryB->domainIndex(elementIndexB);
    m_testTransformations.evaluate(testBasisData, testGeomData, testValues);
  }

  // Iterate over the elements
  for (int indexA = 0; indexA < elementACount; ++indexA) {
    const int elementIndexA = elementIndicesA[indexA];
    rawGeometryA->setupGeometry(elementIndexA, *geometryA);
    if (callVariant == TEST_TRIAL) {
      geometryA->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
      if (testGeomDeps & DOMAIN_INDEX)
        testGeomData.domainIndex = rawGeometryA->domainIndex(elementIndexA);
      m_testTransformations.evaluate(testBasisData, testGeomData, testValues);
    } else {
      geometryA->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
      if (trialGeomDeps & DOMAIN_INDEX)
        trialGeomData.domainIndex = rawGeometryA->domainIndex(elementIndexA);
      m_trialTransformations.evaluate(trialBasisData, trialGeomData,
                                      trialValues);
    }

    m_kernels.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);
    m_integral.evaluateWithNontensorQuadratureRule(
        testGeomData, trialGeomData, testValues, trialValues, kernelValues,
        m_quadWeights, *result[indexA]);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void NonseparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    integrate(const std::vector<ElementIndexPair> &elementIndexPairs,
              const Shapeset<BasisFunctionType> &testShapeset,
              const Shapeset<BasisFunctionType> &trialShapeset,
              const std::vector<arma::Mat<ResultType> *> &result) const {
  const int pointCount = m_quadWeights.size();
  const int geometryPairCount = elementIndexPairs.size();

  if (result.size() != elementIndexPairs.size())
    throw std::invalid_argument(
        "NonseparableNumericalTestKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'elementIndicesA' must have the same number "
        "of elements");
  if (pointCount == 0 || geometryPairCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryPairCount != 0, set elements of result to 0.

  const int testDofCount = testShapeset.size();
  const int trialDofCount = trialShapeset.size();

  const BasisData<BasisFunctionType> &testBasisData =
      basisData(TEST, testShapeset);
  const BasisData<BasisFunctionType> &trialBasisData =
      basisData(TRIAL, trialShapeset);
  GeometricalData<CoordinateType> &testGeomData = m_testGeomData.local();
  GeometricalData<CoordinateType> &trialGeomData = m_trialGeomData.local();

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> testGeometry(m_testGeometryFactory.make());
  std::unique_ptr<Geometry> trialGeometry(m_trialGeometryFactory.make());

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
  CollectionOf3dArrays<KernelType> kernelValues;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->set_size(testDofCount, trialDofCount);
  }

  // Iterate over the elements
  for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex) {
    const int testElementIndex = elementIndexPairs[pairIndex].first;
    const int trialElementIndex = elementIndexPairs[pairIndex].second;
    m_testRawGeometry.setupGeometry(testElementIndex, *testGeometry);
    m_trialRawGeometry.setupGeometry(trialElementIndex, *trialGeometry);
    testGeometry->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
    if (testGeomDeps & DOMAIN_INDEX)
      testGeomData.domainIndex =
          m_testRawGeometry.domainIndex(testElementIndex);
    trialGeometry->getData(trialGeomDeps, m_localTrialQuadPoints,
                           trialGeomData);
    if (trialGeomDeps & DOMAIN_INDEX)
      trialGeomData.domainIndex =
          m_trialRawGeometry.domainIndex(trialElementIndex);
    m_testTransformations.evaluate(testBasisData, testGeomData, testValues);
    m_trialTransformations.evaluate(trialBasisData, trialGeomData, trialValues);

    m_kernels.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);
    m_integral.evaluateWithNontensorQuadratureRule(
        testGeomData, trialGeomData, testValues, trialValues, kernelValues,
        m_quadWeights, *result[pairIndex]);
  }
}

} // namespace Fiber
