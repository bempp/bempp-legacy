#include "separable_numerical_double_integrator.hpp" // To keep IDEs happy

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
#include <memory>

namespace Fiber
{

template <typename ValueType, typename GeometryFactory>
SeparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::
SeparableNumericalDoubleIntegrator(
        const arma::Mat<ValueType>& localTestQuadPoints,
        const arma::Mat<ValueType>& localTrialQuadPoints,
        const std::vector<ValueType> testQuadWeights,
        const std::vector<ValueType> trialQuadWeights,
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
    m_testQuadWeights(testQuadWeights),
    m_trialQuadWeights(trialQuadWeights),
    m_geometryFactory(geometryFactory),
    m_vertices(vertices),
    m_elementCornerIndices(elementCornerIndices),
    m_auxElementData(auxElementData),
    m_testExpression(testExpression),
    m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_openClOptions(openClOptions)
{}

template <typename ValueType, typename GeometryFactory>
void SeparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::setupGeometry(
        int elementIndex, typename GeometryFactory::Geometry& geometry) const
{
    const int dimGrid = m_vertices.n_rows;
    int cornerCount = 0;
    for (; cornerCount < m_elementCornerIndices.n_rows; ++cornerCount)
        if (m_elementCornerIndices(cornerCount, elementIndex) < 0)
            break;
    arma::Mat<ValueType> corners(dimGrid, cornerCount);
    for (int cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex)
        corners.col(cornerIndex) = m_vertices.col(
                    m_elementCornerIndices(cornerIndex, elementIndex));
    geometry.setup(corners, m_auxElementData.unsafe_col(elementIndex));
}

template <typename ValueType, typename GeometryFactory>
void SeparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::integrate(
        CallVariant callVariant,
        const std::vector<int>& elementIndicesA,
        int elementIndexB,
        const Basis<ValueType>& basisA,
        const Basis<ValueType>& basisB,
        LocalDofIndex localDofIndexB,
        arma::Cube<ValueType>& result) const
{
    const int testPointCount = m_localTestQuadPoints.n_cols;
    const int trialPointCount = m_localTrialQuadPoints.n_cols;
    const int elementACount = elementIndicesA.size();

    if (testPointCount == 0 || trialPointCount == 0 || elementACount == 0)
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
    Array4D<ValueType> kernelValues(kernelRowCount, kernelColCount,
                                    testPointCount, trialPointCount);

    result.set_size(testDofCount, trialDofCount, elementACount);

    setupGeometry(elementIndexB, *geometryB);
    arma::Mat<ValueType> corn;
    geometryB->corners(corn);
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
        setupGeometry(elementIndicesA[indexA], *geometryA);
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

        m_kernel.evaluateOnGrid(testGeomData, trialGeomData, kernelValues);

        // For now, we assume that the kernel is (general) tensorial,
        // later we might handle specially the case of it being a scalar
        // times the identity tensor.
        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
            for (int testDof = 0; testDof < testDofCount; ++testDof)
            {
                ValueType sum = 0.;
                for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
                    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
                        for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                            for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                sum +=  m_testQuadWeights[testPoint] *
                                        testGeomData.integrationElements(testPoint) *
                                        testValues(testDim, testDof, testPoint) *
                                        kernelValues(testDim, trialDim, testPoint, trialPoint) *
                                        trialValues(trialDim, trialDof, trialPoint) *
                                        trialGeomData.integrationElements(trialPoint) *
                                        m_trialQuadWeights[trialPoint];
                result(testDof, trialDof, indexA) = sum;
            }
    }
}

template <typename ValueType, typename GeometryFactory>
void SeparableNumericalDoubleIntegrator<ValueType, GeometryFactory>::integrate(
            const std::vector<ElementIndexPair>& elementIndexPairs,
            const Basis<ValueType>& testBasis,
            const Basis<ValueType>& trialBasis,
            arma::Cube<ValueType>& result) const
{
    const int testPointCount = m_localTestQuadPoints.n_cols;
    const int trialPointCount = m_localTrialQuadPoints.n_cols;
    const int geometryPairCount = elementIndexPairs.size();

    if (testPointCount == 0 || trialPointCount == 0 || geometryPairCount == 0)
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
    Array4D<ValueType> kernelValues(kernelRowCount, kernelColCount,
                                    testPointCount, trialPointCount);

    result.set_size(testDofCount, trialDofCount, geometryPairCount);

    testBasis.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);

    // Iterate over the elements
    for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex)
    {
        setupGeometry(elementIndexPairs[pairIndex].first, *testGeometry);
        setupGeometry(elementIndexPairs[pairIndex].first, *trialGeometry);
        testGeometry->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        trialGeometry->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_testExpression.evaluate(testBasisData, testGeomData, testValues);
        m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);

        m_kernel.evaluateOnGrid(testGeomData, trialGeomData, kernelValues);

        // For now, we assume that the kernel is (general) tensorial,
        // later we might handle specially the case of it being a scalar
        // times the identity tensor.
        for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
            for (int testDof = 0; testDof < testDofCount; ++testDof)
            {
                ValueType sum = 0.;
                for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
                    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
                        for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                            for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                sum +=  m_testQuadWeights[testPoint] *
                                        testGeomData.integrationElements(testPoint) *
                                        testValues(testDim, testDof, testPoint) *
                                        kernelValues(testDim, trialDim, testPoint, trialPoint) *
                                        trialValues(trialDim, trialDof, trialPoint) *
                                        trialGeomData.integrationElements(trialPoint) *
                                        m_trialQuadWeights[trialPoint];
                result(testDof, trialDof, pairIndex) = sum;
            }
    }
}


} // namespace Fiber
