#include "nonseparable_numerical_double_integrator.hpp" // To keep IDEs happy

#include "array_2d.hpp"
#include "array_3d.hpp"
#include "array_4d.hpp"

#include "basis.hpp"
#include "basis_data.hpp"
#include "expression.hpp"
#include "geometrical_data.hpp"
#include "kernel.hpp"
#include "opencl_options.hpp"
#include "types.hpp"

#include <cassert>
#include <iostream>
#include <memory>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
NonseparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::
NonseparableNumericalDoubleIntegrator(
        const arma::Mat<ValueType>& localTestQuadPoints,
        const arma::Mat<ValueType>& localTrialQuadPoints,
        const std::vector<ValueType> quadWeights,
        const GeometryFactory& geometryFactory,
        const arma::Mat<ValueType>& vertices,
        const arma::Mat<int>& elementCornerIndices,
        const arma::Mat<char>& auxElementData,
        const Expression<ValueType>& testExpression,
        const Kernel<ValueType>& kernel,
        const Expression<ValueType>& trialExpression,
        const OpenClOptions& openClOptions) :
    m_localTestQuadPoints(localTestQuadPoints),
    m_localTrialQuadPoints(localTrialQuadPoints),
    m_quadWeights(quadWeights),
    m_geometryFactory(geometryFactory),
    m_vertices(vertices),
    m_elementCornerIndices(elementCornerIndices),
    m_auxElementData(auxElementData),
    m_testExpression(testExpression),
    m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_openClOptions(openClOptions)
{
    const int pointCount = quadWeights.size();
    if (localTestQuadPoints.n_cols != pointCount ||
            localTrialQuadPoints.n_cols != pointCount)
        throw std::invalid_argument("NonseparableNumericalDoubleIntegrator::"
                                    "NonseparableNumericalDoubleIntegrator(): "
                                    "numbers of points and weights do not match");
}


template <typename ValueType, typename GeometryFactory>
inline void NonseparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::
setupGeometryConveniently(
        int elementIndex, typename GeometryFactory::Geometry& geometry) const
{
    setupGeometry(elementIndex,
                  m_vertices, m_elementCornerIndices, m_auxElementData,
                  geometry);
}

template <typename ValueType, typename GeometryFactory>
void NonseparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::integrate(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        const Basis<ValueType>& basisA,
        const Basis<ValueType>& basisB,
        LocalDofIndex localDofIndexB,
        arma::Cube<ValueType>& result) const
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

    // TODO: This will need to be modified once we allow scalar-valued kernels
    // (treated as if they were multiplied by the unit tensor) with
    // vector-valued functions
    assert(testComponentCount == kernelRowCount);
    assert(kernelColCount == trialComponentCount);

    BasisData<ValueType> testBasisData, trialBasisData;
    GeometricalData<ValueType> testGeomData, trialGeomData;

    int testBasisDeps = 0, trialBasisDeps = 0;
    int testGeomDeps = INTEGRATION_ELEMENTS;
    int trialGeomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, testGeomDeps);
    m_trialExpression.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> geometryA(m_geometryFactory.make());
    std::auto_ptr<Geometry> geometryB(m_geometryFactory.make());

    arma::Cube<ValueType> testValues, trialValues;
    arma::Cube<ValueType> kernelValues;

    result.set_size(testDofCount, trialDofCount, elementACount);

    setupGeometryConveniently(elementIndexB, *geometryB);
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
        setupGeometryConveniently(elementIndicesA[indexA], *geometryA);
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

        // For now, we assume that the kernel is (general) tensorial,
        // later we might handle specially the case of it being a scalar
        // times the identity tensor.
        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
            for (int testDof = 0; testDof < testDofCount; ++testDof)
            {
                ValueType sum = 0.;
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

template <typename ValueType, typename GeometryFactory>
void NonseparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::integrate(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<ValueType>& testBasis,
            const Basis<ValueType>& trialBasis,
            arma::Cube<ValueType>& result) const
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

    // TODO: This will need to be modified once we allow scalar-valued kernels
    // (treated as if they were multiplied by the unit tensor) with
    // vector-valued functions
    assert(testComponentCount == kernelRowCount);
    assert(kernelColCount == trialComponentCount);

    BasisData<ValueType> testBasisData, trialBasisData;
    GeometricalData<ValueType> testGeomData, trialGeomData;

    int testBasisDeps = 0, trialBasisDeps = 0;
    int testGeomDeps = INTEGRATION_ELEMENTS;
    int trialGeomDeps = INTEGRATION_ELEMENTS;

    m_testExpression.addDependencies(testBasisDeps, testGeomDeps);
    m_trialExpression.addDependencies(trialBasisDeps, trialGeomDeps);
    m_kernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    typedef typename GeometryFactory::Geometry Geometry;
    std::auto_ptr<Geometry> testGeometry(m_geometryFactory.make());
    std::auto_ptr<Geometry> trialGeometry(m_geometryFactory.make());

    arma::Cube<ValueType> testValues, trialValues;
    arma::Cube<ValueType> kernelValues;

    result.set_size(testDofCount, trialDofCount, geometryPairCount);

    testBasis.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex)
    {
        setupGeometryConveniently(elementIndexPairs[pairIndex].first, *testGeometry);
        setupGeometryConveniently(elementIndexPairs[pairIndex].first, *trialGeometry);
        testGeometry->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        trialGeometry->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_testExpression.evaluate(testBasisData, testGeomData, testValues);
        m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);

        m_kernel.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);

        // For now, we assume that the kernel is (general) tensorial,
        // later we might handle specially the case of it being a scalar
        // times the identity tensor.
        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
            for (int testDof = 0; testDof < testDofCount; ++testDof)
            {
                ValueType sum = 0.;
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
