#include "separable_numerical_double_integrator.hpp" // To keep IDEs happy

#include "array_2d.hpp"
#include "array_3d.hpp"
#include "array_4d.hpp"

#include <cassert>

#include "basis.hpp"
#include "basis_data.hpp"
#include "expression.hpp"
#include "geometrical_data.hpp"
#include "kernel.hpp"
#include "opencl_options.hpp"
#include "types.hpp"

namespace Fiber
{

template <typename ValueType, typename GeometryImp>
SeparableNumericalDoubleIntegrator<ValueType, GeometryImp>::
SeparableNumericalDoubleIntegrator(
        const arma::Mat<ValueType>& localTestQuadPoints,
        const arma::Mat<ValueType>& localTrialQuadPoints,
        const std::vector<ValueType> testQuadWeights,
        const std::vector<ValueType> trialQuadWeights,
        const Expression<ValueType>& testExpression,
        const Kernel<ValueType>& kernel,
        const Expression<ValueType>& trialExpression,
        const OpenClOptions& openClOptions) :
    m_localTestQuadPoints(localTestQuadPoints),
    m_localTrialQuadPoints(localTrialQuadPoints),
    m_testQuadWeights(testQuadWeights),
    m_trialQuadWeights(trialQuadWeights),
    m_testExpression(testExpression),
    m_kernel(kernel),
    m_trialExpression(trialExpression),
    m_openClOptions(openClOptions)
{}

template <typename ValueType, typename GeometryImp>
void SeparableNumericalDoubleIntegrator<ValueType, GeometryImp>::integrate(
        CallVariant callVariant,
        const std::vector<const GeometryImp*>& geometriesA,
        const GeometryImp& geometryB,
        const Basis<ValueType>& basisA,
        const Basis<ValueType>& basisB,
        LocalDofIndex localDofIndexB,
        arma::Cube<ValueType>& result) const
{
    const int testPointCount = m_localTestQuadPoints.n_cols;
    const int trialPointCount = m_localTrialQuadPoints.n_cols;
    const int geometryACount = geometriesA.size();

    if (testPointCount == 0 || trialPointCount == 0 || geometryACount == 0)
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

    arma::Cube<ValueType> testValues, trialValues;
    Array4D<ValueType> kernelValues(kernelRowCount, kernelColCount,
                                    testPointCount, trialPointCount);

    result.set_size(testDofCount, trialDofCount, geometryACount);

    if (callVariant == TEST_TRIAL)
    {
        basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
        basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB, trialBasisData);
        geometryB.getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
        m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);
    }
    else
    {
        basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);
        basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB, testBasisData);
        geometryB.getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        m_testExpression.evaluate(testBasisData, testGeomData, testValues);
    }

    // Iterate over the elements
    for (int indexA = 0; indexA < geometryACount; ++indexA)
    {
        const GeometryImp& geometryA = *geometriesA[indexA];
        if (callVariant == TEST_TRIAL)
        {
            geometryA.getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
            m_testExpression.evaluate(testBasisData, testGeomData, testValues);
        }
        else
        {
            geometryA.getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
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

template <typename ValueType, typename GeometryImp>
void SeparableNumericalDoubleIntegrator<ValueType, GeometryImp>::integrate(
            const std::vector<GeometryImpPair>& geometryPairs,
            const Basis<ValueType>& testBasis,
            const Basis<ValueType>& trialBasis,
            arma::Cube<ValueType>& result) const
{
    const int testPointCount = m_localTestQuadPoints.n_cols;
    const int trialPointCount = m_localTrialQuadPoints.n_cols;
    const int geometryPairCount = geometryPairs.size();

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

    arma::Cube<ValueType> testValues, trialValues;
    Array4D<ValueType> kernelValues(kernelRowCount, kernelColCount,
                                    testPointCount, trialPointCount);

    result.set_size(testDofCount, trialDofCount, geometryPairCount);

    testBasis.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS, testBasisData);
    trialBasis.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS, trialBasisData);
//        geometryB.getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
    //    m_trialExpression.evaluate(trialBasisData, trialGeomData, trialValues);

    // Iterate over the elements
    for (int pairIndex = 0; pairIndex < geometryPairCount; ++pairIndex)
    {
        const GeometryImp& testGeometry = *geometryPairs[pairIndex].first;
        const GeometryImp& trialGeometry = *geometryPairs[pairIndex].second;
        testGeometry.getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
        m_testExpression.evaluate(testBasisData, testGeomData, testValues);
        trialGeometry.getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
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
