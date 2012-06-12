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

#include "array_2d.hpp"
#include "array_3d.hpp"
#include "array_4d.hpp"

#include "basis.hpp"
#include "basis_data.hpp"
#include "conjugate.hpp"
#include "expression.hpp"
#include "geometrical_data.hpp"
#include "kernel.hpp"
#include "opencl_handler.hpp"
#include "raw_grid_geometry.hpp"
#include "types.hpp"

#include <cassert>
#include <iostream>
#include <memory>

namespace Fiber
{

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
NonseparableNumericalTestKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
NonseparableNumericalTestKernelTrialIntegrator(
        const arma::Mat<CoordinateType>& localTestQuadPoints,
        const arma::Mat<CoordinateType>& localTrialQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const GeometryFactory& testGeometryFactory,
        const GeometryFactory& trialGeometryFactory,
        const RawGridGeometry<CoordinateType>& testRawGeometry,
        const RawGridGeometry<CoordinateType>& trialRawGeometry,
        const ExpressionList<ResultType>& testExpressionList,
        const Kernel<KernelType>& kernel,
        const ExpressionList<ResultType>& trialExpressionList,
        const OpenClHandler& openClHandler) :
    m_localTestQuadPoints(localTestQuadPoints),
    m_localTrialQuadPoints(localTrialQuadPoints),
    m_quadWeights(quadWeights),
    m_testGeometryFactory(testGeometryFactory),
    m_trialGeometryFactory(trialGeometryFactory),
    m_testRawGeometry(testRawGeometry),
    m_trialRawGeometry(trialRawGeometry),
    m_testExpressionList(testExpressionList),
    m_kernel(kernel),
    m_trialExpressionList(trialExpressionList),
    m_openClHandler(openClHandler)
{
    const size_t pointCount = quadWeights.size();
    if (localTestQuadPoints.n_cols != pointCount ||
            localTrialQuadPoints.n_cols != pointCount)
        throw std::invalid_argument("NonseparableNumericalTestKernelTrialIntegrator::"
                                    "NonseparableNumericalTestKernelTrialIntegrator(): "
                                    "numbers of points and weights do not match");

    const size_t expressionCount = testExpressionList.termCount();
    if (expressionCount != trialExpressionList.termCount())
        throw std::invalid_argument("SeparableNumericalTestKernelTrialIntegrator::"
                                    "SeparableNumericalTestKernelTrialIntegrator(): "
                                    "test and trial expression lists have "
                                    "different lengths");
    assert(expressionCount > 0);
    // Multiply the test and trial expression weigths and store them
    m_expressionWeights.resize(expressionCount);
    for (size_t i = 0; i < expressionCount; ++i)
        m_expressionWeights[i] =
                testExpressionList.weight(i) * trialExpressionList.weight(i);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
NonseparableNumericalTestKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
integrate(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        const Basis<BasisFunctionType>& basisA,
        const Basis<BasisFunctionType>& basisB,
        LocalDofIndex localDofIndexB,
        arma::Cube<ResultType>& result) const
{
    const int pointCount = m_quadWeights.size();
    const int elementACount = elementIndicesA.size();

    if (pointCount == 0 || elementACount == 0)
        return;
    // TODO: in the (pathological) case that pointCount == 0 but
    // geometryCount != 0, set elements of result to 0.

    // Evaluate constants
    const int testComponentCount = m_testExpressionList.codomainDimension();
    const int trialComponentCount = m_trialExpressionList.codomainDimension();
    // In the constructor we have checked that m_trialExpressionList has
    // the same number of terms
    const int expressionCount = m_testExpressionList.termCount();

    const int dofCountA = basisA.size();
    const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
    const int testDofCount = callVariant == TEST_TRIAL ? dofCountA : dofCountB;
    const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB : dofCountA;

    const int kernelRowCount = m_kernel.codomainDimension();
    const int kernelColCount = m_kernel.domainDimension();

    // Assert that the kernel tensor dimensions are compatible
    // with the number of components of the functions

    const bool scalarKernel = (kernelRowCount == 1 && kernelColCount == 1);
    if (scalarKernel)
        assert(testComponentCount == trialComponentCount);
    else
    {
        assert(testComponentCount == kernelRowCount);
        assert(kernelColCount == trialComponentCount);
    }

    BasisData<BasisFunctionType> testBasisData, trialBasisData;
    GeometricalData<CoordinateType> testGeomData, trialGeomData;

    size_t testBasisDeps = 0, trialBasisDeps = 0;
    size_t testGeomDeps = INTEGRATION_ELEMENTS;
    size_t trialGeomDeps = INTEGRATION_ELEMENTS;

    m_testExpressionList.addDependencies(testBasisDeps, testGeomDeps);
    m_trialExpressionList.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;

    std::auto_ptr<Geometry> geometryA, geometryB;
    const RawGridGeometry<CoordinateType> *rawGeometryA = 0, *rawGeometryB = 0;
    if (callVariant == TEST_TRIAL)
    {
        geometryA = m_testGeometryFactory.make();
        geometryB = m_trialGeometryFactory.make();
        rawGeometryA = &m_testRawGeometry;
        rawGeometryB = &m_trialRawGeometry;
    }
    else
    {
        geometryA = m_trialGeometryFactory.make();
        geometryB = m_testGeometryFactory.make();
        rawGeometryA = &m_trialRawGeometry;
        rawGeometryB = &m_testRawGeometry;
    }

    std::vector<arma::Cube<BasisFunctionType> > testValues, trialValues;
    arma::Cube<KernelType> kernelValues;

    result.set_size(testDofCount, trialDofCount, elementACount);

    rawGeometryB->setupGeometry(elementIndexB, *geometryB);
    if (callVariant == TEST_TRIAL)
    {
        basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
        basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB, trialBasisData);
        geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_trialExpressionList.evaluate(trialBasisData, trialGeomData, trialValues);
    }
    else
    {
        basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);
        basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB, testBasisData);
        geometryB->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        m_testExpressionList.evaluate(testBasisData, testGeomData, testValues);
    }

    // Iterate over the elements
    for (int indexA = 0; indexA < elementACount; ++indexA)
    {
        rawGeometryA->setupGeometry(elementIndicesA[indexA], *geometryA);
        if (callVariant == TEST_TRIAL)
        {
            geometryA->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
            m_testExpressionList.evaluate(testBasisData, testGeomData, testValues);
        }
        else
        {
            geometryA->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
            m_trialExpressionList.evaluate(trialBasisData, trialGeomData, trialValues);
        }

        m_kernel.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);

        if (scalarKernel)
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int expression = 0; expression < expressionCount; ++expression)
                            for (int dim = 0; dim < testComponentCount; ++dim)
                                sum +=  m_quadWeights[point] *
                                        testGeomData.integrationElements(point) *
                                        m_expressionWeights[expression] *
                                        conjugate(testValues[expression](dim, testDof, point)) *
                                        kernelValues(0, 0, point) *
                                        trialValues[expression](dim, trialDof, point) *
                                        trialGeomData.integrationElements(point);
                    result(testDof, trialDof, indexA) = sum;
                }
        else
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int expression = 0; expression < expressionCount; ++expression)
                            for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                                for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                    sum +=  m_quadWeights[point] *
                                            testGeomData.integrationElements(point) *
                                            m_expressionWeights[expression] *
                                            conjugate(testValues[expression](testDim, testDof, point)) *
                                            kernelValues(testDim, trialDim, point) *
                                            trialValues[expression](trialDim, trialDof, point) *
                                            trialGeomData.integrationElements(point);
                    result(testDof, trialDof, indexA) = sum;
                }
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
NonseparableNumericalTestKernelTrialIntegrator<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::
integrate(
        const std::vector<ElementIndexPair>& elementIndexPairs,
        const Basis<BasisFunctionType>& testBasis,
        const Basis<BasisFunctionType>& trialBasis,
        arma::Cube<ResultType>& result) const
{
    const int pointCount = m_quadWeights.size();
    const int geometryPairCount = elementIndexPairs.size();

    if (pointCount == 0 || geometryPairCount == 0)
        return;
    // TODO: in the (pathological) case that pointCount == 0 but
    // geometryPairCount != 0, set elements of result to 0.

    // Evaluate constants
    const int testComponentCount = m_testExpressionList.codomainDimension();
    const int trialComponentCount = m_trialExpressionList.codomainDimension();
    // In the constructor we have checked that m_trialExpressionList has
    // the same number of terms
    const int expressionCount = m_testExpressionList.termCount();

    const int testDofCount = testBasis.size();
    const int trialDofCount = trialBasis.size();

    const int kernelRowCount = m_kernel.codomainDimension();
    const int kernelColCount = m_kernel.domainDimension();

    // Assert that the kernel tensor dimensions are compatible
    // with the number of components of the functions

    const bool scalarKernel = (kernelRowCount == 1 && kernelColCount == 1);
    if (scalarKernel)
        assert(testComponentCount == trialComponentCount);
    else
    {
        assert(testComponentCount == kernelRowCount);
        assert(kernelColCount == trialComponentCount);
    }

    BasisData<BasisFunctionType> testBasisData, trialBasisData;
    GeometricalData<CoordinateType> testGeomData, trialGeomData;

    size_t testBasisDeps = 0, trialBasisDeps = 0;
    size_t testGeomDeps = INTEGRATION_ELEMENTS;
    size_t trialGeomDeps = INTEGRATION_ELEMENTS;

    m_testExpressionList.addDependencies(testBasisDeps, testGeomDeps);
    m_trialExpressionList.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> testGeometry(m_testGeometryFactory.make());
    std::auto_ptr<Geometry> trialGeometry(m_trialGeometryFactory.make());

    std::vector<arma::Cube<BasisFunctionType> > testValues, trialValues;
    arma::Cube<KernelType> kernelValues;

    result.set_size(testDofCount, trialDofCount, geometryPairCount);

    testBasis.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex)
    {
        m_testRawGeometry.setupGeometry(elementIndexPairs[pairIndex].first, *testGeometry);
        m_trialRawGeometry.setupGeometry(elementIndexPairs[pairIndex].second, *trialGeometry);
        testGeometry->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        trialGeometry->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_testExpressionList.evaluate(testBasisData, testGeomData, testValues);
        m_trialExpressionList.evaluate(trialBasisData, trialGeomData, trialValues);

        m_kernel.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);

        if (scalarKernel)
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int expression = 0; expression < expressionCount; ++expression)
                            for (int dim = 0; dim < testComponentCount; ++dim)
                                sum +=  m_quadWeights[point] *
                                        testGeomData.integrationElements(point) *
                                        m_expressionWeights[expression] *
                                        conjugate(testValues[expression](dim, testDof, point)) *
                                        kernelValues(0, 0, point) *
                                        trialValues[expression](dim, trialDof, point) *
                                        trialGeomData.integrationElements(point);
                    result(testDof, trialDof, pairIndex) = sum;
                }
        else
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int expression = 0; expression < expressionCount; ++expression)
                            for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                                for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                    sum +=  m_quadWeights[point] *
                                            testGeomData.integrationElements(point) *
                                            m_expressionWeights[expression] *
                                            conjugate(testValues[expression](testDim, testDof, point)) *
                                            kernelValues(testDim, trialDim, point) *
                                            trialValues[expression](trialDim, trialDof, point) *
                                            trialGeomData.integrationElements(point);
                    result(testDof, trialDof, pairIndex) = sum;
                }

    }
}

} // namespace Fiber
