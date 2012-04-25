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

#include "nonseparable_numerical_test_kernel_trial_integrator.hpp" // To keep IDEs happy

#include "array_2d.hpp"
#include "array_3d.hpp"
#include "array_4d.hpp"

#include "basis.hpp"
#include "basis_data.hpp"
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

template <typename BasisValueType, typename KernelValueType, typename GeometryFactory>
NonseparableNumericalTestKernelTrialIntegrator<BasisValueType, KernelValueType, GeometryFactory>::
NonseparableNumericalTestKernelTrialIntegrator(
        const arma::Mat<CoordinateType>& localTestQuadPoints,
        const arma::Mat<CoordinateType>& localTrialQuadPoints,
        const std::vector<CoordinateType> quadWeights,
        const GeometryFactory& geometryFactory,
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const Expression<BasisValueType>& testExpression,
        const Kernel<KernelValueType>& kernel,
        const Expression<BasisValueType>& trialExpression,
        const OpenClHandler<CoordinateType, int>& openClHandler) :
    m_localTestQuadPoints(localTestQuadPoints),
    m_localTrialQuadPoints(localTrialQuadPoints),
    m_quadWeights(quadWeights),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_testExpression(testExpression),
    m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_openClHandler(openClHandler)
{
    const int pointCount = quadWeights.size();
    if (localTestQuadPoints.n_cols != pointCount ||
            localTrialQuadPoints.n_cols != pointCount)
        throw std::invalid_argument("NonseparableNumericalTestKernelTrialIntegrator::"
                                    "NonseparableNumericalTestKernelTrialIntegrator(): "
                                    "numbers of points and weights do not match");
}

template <typename BasisValueType, typename KernelValueType, typename GeometryFactory>
void NonseparableNumericalTestKernelTrialIntegrator<BasisValueType, KernelValueType, GeometryFactory>::integrate(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        const Basis<BasisValueType>& basisA,
        const Basis<BasisValueType>& basisB,
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
    const int testComponentCount = m_testExpression.codomainDimension();
    const int trialComponentCount = m_trialExpression.codomainDimension();
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

    BasisData<BasisValueType> testBasisData, trialBasisData;
    GeometricalData<CoordinateType> testGeomData, trialGeomData;

    int testBasisDeps = 0, trialBasisDeps = 0;
    int testGeomDeps = INTEGRATION_ELEMENTS;
    int trialGeomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, testGeomDeps);
    m_trialExpression.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometryA(m_geometryFactory.make());
    std::auto_ptr<Geometry> geometryB(m_geometryFactory.make());

    arma::Cube<BasisValueType> testValues, trialValues;
    arma::Cube<KernelValueType> kernelValues;

    result.set_size(testDofCount, trialDofCount, elementACount);

    m_rawGeometry.setupGeometry(elementIndexB, *geometryB);
    if (callVariant == TEST_TRIAL)
    {
        basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
        basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB, trialBasisData);
        geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);
    }
    else
    {
        basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);
        basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB, testBasisData);
        geometryB->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        m_testExpression.evaluate(testBasisData, testGeomData, testValues);
    }

    // Iterate over the elements
    for (int indexA = 0; indexA < elementACount; ++indexA)
    {
        m_rawGeometry.setupGeometry(elementIndicesA[indexA], *geometryA);
        if (callVariant == TEST_TRIAL)
        {
            geometryA->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
            m_testExpression.evaluate(testBasisData, testGeomData, testValues);
        }
        else
        {
            geometryA->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
            m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);
        }

        m_kernel.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);

        if (scalarKernel)
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int dim = 0; dim < testComponentCount; ++dim)
                            sum +=  m_quadWeights[point] *
                                    testGeomData.integrationElements(point) *
                                    testValues(dim, testDof, point) *
                                    kernelValues(0, 0, point) *
                                    trialValues(dim, trialDof, point) *
                                    trialGeomData.integrationElements(point);
                    result(testDof, trialDof, indexA) = sum;
                }
        else
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                            for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                sum +=  m_quadWeights[point] *
                                        testGeomData.integrationElements(point) *
                                        testValues(testDim, testDof, point) *
                                        kernelValues(testDim, trialDim, point) *
                                        trialValues(trialDim, trialDof, point) *
                                        trialGeomData.integrationElements(point);
                    result(testDof, trialDof, indexA) = sum;
                }
    }
}

template <typename BasisValueType, typename KernelValueType, typename GeometryFactory>
void NonseparableNumericalTestKernelTrialIntegrator<BasisValueType, KernelValueType, GeometryFactory>::integrate(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<BasisValueType>& testBasis,
            const Basis<BasisValueType>& trialBasis,
            arma::Cube<ResultType>& result) const
{
    const int pointCount = m_quadWeights.size();
    const int geometryPairCount = elementIndexPairs.size();

    if (pointCount == 0 || geometryPairCount == 0)
        return;
    // TODO: in the (pathological) case that pointCount == 0 but
    // geometryPairCount != 0, set elements of result to 0.

    // Evaluate constants
    const int testComponentCount = m_testExpression.codomainDimension();
    const int trialComponentCount = m_trialExpression.codomainDimension();
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

    BasisData<BasisValueType> testBasisData, trialBasisData;
    GeometricalData<CoordinateType> testGeomData, trialGeomData;

    int testBasisDeps = 0, trialBasisDeps = 0;
    int testGeomDeps = INTEGRATION_ELEMENTS;
    int trialGeomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, testGeomDeps);
    m_trialExpression.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> testGeometry(m_geometryFactory.make());
    std::auto_ptr<Geometry> trialGeometry(m_geometryFactory.make());

    arma::Cube<BasisValueType> testValues, trialValues;
    arma::Cube<KernelValueType> kernelValues;

    result.set_size(testDofCount, trialDofCount, geometryPairCount);

    testBasis.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex)
    {
        m_rawGeometry.setupGeometry(elementIndexPairs[pairIndex].first, *testGeometry);
        m_rawGeometry.setupGeometry(elementIndexPairs[pairIndex].second, *trialGeometry);
        testGeometry->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        trialGeometry->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_testExpression.evaluate(testBasisData, testGeomData, testValues);
        m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);

        m_kernel.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);

        if (scalarKernel)
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int dim = 0; dim < testComponentCount; ++dim)
                            sum +=  m_quadWeights[point] *
                                    testGeomData.integrationElements(point) *
                                    testValues(dim, testDof, point) *
                                    kernelValues(0, 0, point) *
                                    trialValues(dim, trialDof, point) *
                                    trialGeomData.integrationElements(point);
                    result(testDof, trialDof, pairIndex) = sum;
                }
        else
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ResultType sum = 0.;
                    for (int point = 0; point < pointCount; ++point)
                        for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                            for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                sum +=  m_quadWeights[point] *
                                        testGeomData.integrationElements(point) *
                                        testValues(testDim, testDof, point) *
                                        kernelValues(testDim, trialDim, point) *
                                        trialValues(trialDim, trialDof, point) *
                                        trialGeomData.integrationElements(point);
                    result(testDof, trialDof, pairIndex) = sum;
                }

    }
}

} // namespace Fiber
