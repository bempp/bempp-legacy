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

#include "cuda_separable_numerical_test_kernel_trial_integrator.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
CudaSeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                ResultType, GeometryFactory>::
    CudaSeparableNumericalTestKernelTrialIntegrator(
        const Matrix<CoordinateType> &localTestQuadPoints,
        const Matrix<CoordinateType> &localTrialQuadPoints,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testQuadWeights(testQuadWeights), m_trialQuadWeights(trialQuadWeights) {

  if (localTestQuadPoints.cols() != testQuadWeights.size())
    throw std::invalid_argument(
        "CudaSeparableNumericalTestKernelTrialIntegrator::"
        "CudaSeparableNumericalTestKernelTrialIntegrator(): "
        "numbers of test points and weights do not match");

  if (localTrialQuadPoints.cols() != trialQuadWeights.size())
    throw std::invalid_argument(
        "CudaSeparableNumericalTestKernelTrialIntegrator::"
        "CudaSeparableNumericalTestKernelTrialIntegrator(): "
        "numbers of trial points and weights do not match");
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
CudaSeparableNumericalTestKernelTrialIntegrator<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::~CudaSeparableNumericalTestKernelTrialIntegrator() { }

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void CudaSeparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                                     ResultType, GeometryFactory>::
    integrate(CallVariant callVariant, const std::vector<int> &elementIndicesA,
              int elementIndexB, const Shapeset<BasisFunctionType> &basisA,
              const Shapeset<BasisFunctionType> &basisB,
              LocalDofIndex localDofIndexB,
              const std::vector<Matrix<ResultType> *> &result) const {

//  const int testPointCount = m_localTestQuadPoints.cols();
//  const int trialPointCount = m_localTrialQuadPoints.cols();
//  const int elementACount = elementIndicesA.size();
//
//  if (result.size() != elementIndicesA.size())
//    throw std::invalid_argument(
//        "SeparableNumericalTestKernelTrialIntegrator::integrate(): "
//        "arrays 'result' and 'elementIndicesA' must have the same number "
//        "of elements");
//  if (testPointCount == 0 || trialPointCount == 0 || elementACount == 0)
//    return;
//  // TODO: in the (pathological) case that pointCount == 0 but
//  // geometryCount != 0, set elements of result to 0.
//
//  // Evaluate constants
//
//  const int dofCountA = basisA.size();
//  const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
//  const int testDofCount = callVariant == TEST_TRIAL ? dofCountA : dofCountB;
//  const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB : dofCountA;
//
//  BasisData<BasisFunctionType> testBasisData, trialBasisData;
//  GeometricalData<CoordinateType> *testGeomData = &m_testGeomData.local();
//  GeometricalData<CoordinateType> *trialGeomData = &m_trialGeomData.local();
//  const GeometricalData<CoordinateType> *constTestGeomData = testGeomData;
//  const GeometricalData<CoordinateType> *constTrialGeomData = trialGeomData;
//
//  size_t testBasisDeps = 0, trialBasisDeps = 0;
//  size_t testGeomDeps = 0, trialGeomDeps = 0;
//
//  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
//  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
//  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
//  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
//
//  typedef typename GeometryFactory::Geometry Geometry;
//  std::unique_ptr<Geometry> geometryA, geometryB;
//  const RawGridGeometry<CoordinateType> *rawGeometryA = 0, *rawGeometryB = 0;
//  if (!m_cacheGeometricalData) {
//    if (callVariant == TEST_TRIAL) {
//      geometryA = m_testGeometryFactory.make();
//      geometryB = m_trialGeometryFactory.make();
//      rawGeometryA = &m_testRawGeometry;
//      rawGeometryB = &m_trialRawGeometry;
//    } else {
//      geometryA = m_trialGeometryFactory.make();
//      geometryB = m_testGeometryFactory.make();
//      rawGeometryA = &m_trialRawGeometry;
//      rawGeometryB = &m_testRawGeometry;
//    }
//  }
//
//  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
//  CollectionOf4dArrays<KernelType> kernelValues;
//
//  for (size_t i = 0; i < result.size(); ++i) {
//    assert(result[i]);
//    result[i]->resize(testDofCount, trialDofCount);
//  }
//
//  if (!m_cacheGeometricalData)
//    rawGeometryB->setupGeometry(elementIndexB, *geometryB);
//  if (callVariant == TEST_TRIAL) {
//    basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
//                    testBasisData);
//    basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB,
//                    trialBasisData);
//    if (m_cacheGeometricalData)
//      constTrialGeomData = &m_cachedTrialGeomData[elementIndexB];
//    else {
//      geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, *trialGeomData);
//      if (trialGeomDeps & DOMAIN_INDEX)
//        trialGeomData->domainIndex = rawGeometryB->domainIndex(elementIndexB);
//    }
//    m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
//                                    trialValues);
//  } else {
//    basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
//                    trialBasisData);
//    basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB,
//                    testBasisData);
//    if (m_cacheGeometricalData)
//      constTestGeomData = &m_cachedTestGeomData[elementIndexB];
//    else {
//      geometryB->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
//      if (testGeomDeps & DOMAIN_INDEX)
//        testGeomData->domainIndex = rawGeometryB->domainIndex(elementIndexB);
//    }
//    m_testTransformations.evaluate(testBasisData, *constTestGeomData,
//                                   testValues);
//  }
//
//  // Iterate over the elements
//  for (int indexA = 0; indexA < elementACount; ++indexA) {
//    const int elementIndexA = elementIndicesA[indexA];
//    if (!m_cacheGeometricalData)
//      rawGeometryA->setupGeometry(elementIndexA, *geometryA);
//    if (callVariant == TEST_TRIAL) {
//      if (m_cacheGeometricalData)
//        constTestGeomData = &m_cachedTestGeomData[elementIndicesA[indexA]];
//      else {
//        geometryA->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
//        if (testGeomDeps & DOMAIN_INDEX)
//          testGeomData->domainIndex = rawGeometryA->domainIndex(elementIndexA);
//      }
//      m_testTransformations.evaluate(testBasisData, *constTestGeomData,
//                                     testValues);
//    } else {
//      if (m_cacheGeometricalData)
//        constTrialGeomData = &m_cachedTrialGeomData[elementIndicesA[indexA]];
//      else {
//        geometryA->getData(trialGeomDeps, m_localTrialQuadPoints,
//                           *trialGeomData);
//        if (trialGeomDeps & DOMAIN_INDEX)
//          trialGeomData->domainIndex = rawGeometryA->domainIndex(elementIndexA);
//      }
//      m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
//                                      trialValues);
//    }
//
//    m_kernels.evaluateOnGrid(*constTestGeomData, *constTrialGeomData,
//                             kernelValues);
//    m_integral.evaluateWithTensorQuadratureRule(
//        *constTestGeomData, *constTrialGeomData, testValues, trialValues,
//        kernelValues, m_testQuadWeights, m_trialQuadWeights, *result[indexA]);
//  }
}

} // namespace Fiber
