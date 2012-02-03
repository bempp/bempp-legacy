#include "double_integrator.hpp"

#include "array_2d.hpp"
#include "array_3d.hpp"
#include "array_4d.hpp"

#include <cassert>

namespace Fiber
{

template <typename ValueType,
          typename CoordinateType,
          typename GeometryImp,
          typename TestFunctionFamilyImp,
          typename TrialFunctionFamilyImp,
          typename KernelImp>
void DoubleIntegrator::integrate(
        int pointCount,
        const CoordinateType* localTestQuadPoints,
        const CoordinateType* localTrialQuadPoints,
        const ValueType* quadWeights,
        int testGeometryCount,
        const GeometryImp** testGeometries,
        int trialGeometryCount,
        const GeometryImp** trialGeometries,
        TestFunctionFamilyImp& testFamily,
        TrialFunctionFamilyImp& trialFamily,
        const KernelImp& kernel,
        ValueType* result)
{
    if (pointCount == 0 || testGeometryCount == 0 || trialGeometryCount == 0)
        return;

    // Evaluate constants
    const int elementDim = trialGeometries[0]->dimension();
    const int worldDim = trialGeometries[0]->worldDimension();

    const int testComponentCount = testFamily.codomainDimension();
    const int trialComponentCount = trialFamily.codomainDimension();
    const int testDofCount = testFamily.size();
    const int trialDofCount = trialFamily.size();

    const int kernelRowCount = kernel.codomainDimension();
    const int kernelColCount = kernel.domainDimension();

    // Assert that the kernel tensor dimensions are compatible
    // with the number of components of the functions

    // TODO: This will need to be modified once we allow scalar-valued kernels
    // (treated as if they were multiplied by the unit tensor) with
    // vector-valued functions
    assert(testComponentCount == kernelRowCount);
    assert(kernelColCount == trialComponentCount);

    // Declare containers
    Array2D<CoordinateType> globalTestQuadPoints (worldDim, pointCount);
    Array2D<CoordinateType> globalTrialQuadPoints(worldDim, pointCount);

    Array2D<CoordinateType> testNormals;
    Array2D<CoordinateType> trialNormals;
    if (kernel.needsTestNormal())
        testNormals.set_size(worldDim, pointCount);
    if (kernel.needsTrialNormal())
        trialNormals.set_size(worldDim, pointCount);

    std::vector<CoordinateType> testIntegrationElements (pointCount);
    std::vector<CoordinateType> trialIntegrationElements(pointCount);

    Array3D<CoordinateType> testJT;
    Array3D<CoordinateType> trialJT;
//    if (testFamily.needsJacobianInverseTransposed())
//        testJT.set_size(elementDim, worldDim, pointCount);
//    if (trialFamily.needsJacobianInverseTransposed())
//        trialJT.set_size(elementDim, worldDim, pointCount);

    Array3D<CoordinateType> testJInvT;
    Array3D<CoordinateType> trialJInvT;
    if (testFamily.needsJacobianInverseTransposed())
        testJInvT.set_size(worldDim, elementDim, pointCount);
    if (trialFamily.needsJacobianInverseTransposed())
        trialJInvT.set_size(worldDim, elementDim, pointCount);

    Array3D<ValueType> testFunctionValues(testComponentCount,
                                          testDofCount, pointCount);
    Array3D<ValueType> trialFunctionValues(trialComponentCount,
                                           trialDofCount, pointCount);

    Array3D<ValueType> kernelValues(kernelRowCount, kernelColCount, pointCount);

    Array4D<ValueType> resultArray(testDofCount, testGeometryCount,
                                   trialDofCount, trialGeometryCount, result);

    // Initialize the function families
    testFamily.setEvaluationPoints(pointCount, localTestQuadPoints);
    trialFamily.setEvaluationPoints(pointCount, localTrialQuadPoints);

    // Iterate over the trial elements
    for (int trialIndex = 0; trialIndex < trialGeometryCount; ++trialIndex)
    {
        const GeometryImp* trialGeometry = trialGeometries[trialIndex];

        // Gather data related to the trial element.
        trialGeometry->getInformation(pointCount,
                                      localTrialQuadPoints,
                                      globalTrialQuadPoints.begin(),
                                      trialNormals.begin(),
                                      &*trialIntegrationElements.begin(),
                                      trialJT.begin(),
                                      trialJInvT.begin());
        trialFamily.evaluate(trialJInvT.begin(), trialFunctionValues.begin());

        // Iterate over the test elements
        for (int testIndex = 0; testIndex < testGeometryCount; ++testIndex)
        {
            const GeometryImp* testGeometry = testGeometries[testIndex];

            // Gather data related to the test element.
            testGeometry->getInformation(pointCount,
                                         localTestQuadPoints,
                                         globalTestQuadPoints.begin(),
                                         testNormals.begin(),
                                         &*testIntegrationElements.begin(),
                                         testJT.begin(),
                                         testJInvT.begin());
            testFamily.evaluate(testJInvT.begin(), testFunctionValues.begin());

            // Evaluate kernel
            kernel.evaluate(pointCount,
                            globalTestQuadPoints.begin(),
                            globalTrialQuadPoints.begin(),
                            testNormals.begin(),
                            trialNormals.begin(),
                            kernelValues.begin());

            // For now, we assume that the kernel is (general) tensorial,
            // later we might handle specially the case of it being a scalar
            // times the identity tensor.
            for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                for (int testDof = 0; testDof < testDofCount; ++testDof)
                {
                    ValueType sum = 0.;
                    for (int quadPoint = 0; quadPoint < pointCount; ++quadPoint)
                        for (int trialDim = 0; trialDim < trialComponentCount; ++trialDim)
                            for (int testDim = 0; testDim < testComponentCount; ++testDim)
                                sum +=  testIntegrationElements[quadPoint] *
                                        testFunctionValues(testDim, testDof, quadPoint) *
                                        kernelValues(testDim, trialDim, quadPoint) *
                                        trialFunctionValues(trialDim, trialDof, quadPoint) *
                                        trialIntegrationElements[quadPoint] *
                                        quadWeights[quadPoint];
                    resultArray(testDof, testIndex, trialDof, trialIndex) = sum;
                }
        }
    }
}

} // namespace Fiber
