#include "standard_test_kernel_trial_integral.hpp"

#include "geometrical_data.hpp"

#include <cassert>

namespace Fiber
{

template <typename IntegrandFunctor>
void StandardTestKernelTrialIntegral<IntegrandFunctor>::
addGeometricalDependencies(int& testGeomDeps, int& trialGeomDeps) const
{
    testGeomDeps |= INTEGRATION_ELEMENTS;
    trialGeomDeps |= INTEGRATION_ELEMENTS;

    m_functor.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
}

template <typename IntegrandFunctor>
void StandardTestKernelTrialIntegral<IntegrandFunctor>::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType>& testValues,
        const CollectionOf3dArrays<BasisFunctionType>& trialValues,
        const CollectionOf4dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& testQuadWeights,
        const std::vector<CoordinateType>& trialQuadWeights,
        arma::Mat<ResultType>& result) const
{
    // Evaluate constants

    const int testDofCount = testValues[0].extent(1);
    const int trialDofCount = trialValues[0].extent(1);

    const int testPointCount = testQuadWeights.size();
    const int trialPointCount = trialQuadWeights.size();

    // Assert that array dimensions are correct

    for (int i = 0; i < kernelValues.size(); ++i) {
        assert(kernelValues[i].extent(2) == testPointCount);
        assert(kernelValues[i].extent(3) == trialPointCount);
    }
    for (int i = 0; i < testValues.size(); ++i)
        assert(testValues[i].extent(2) == testPointCount);
    for (int i = 0; i < trialValues.size(); ++i)
        assert(trialValues[i].extent(2) == trialPointCount);

    // Integrate

    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (int testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
                for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
                    sum += m_functor.evaluate(
                                testGeomData.const_slice(testPoint),
                                trialGeomData.const_slice(trialPoint),
                                testValues.const_slice(testDof, testPoint),
                                trialValues.const_slice(trialDof, trialPoint),
                                kernelValues.const_slice(testPoint, trialPoint)) *
                            testGeomData.integrationElements(testPoint) *
                            trialGeomData.integrationElements(trialPoint) *
                            testQuadWeights[testPoint] *
                            trialQuadWeights[trialPoint];
            result(testDof, trialDof) = sum;
        }
}

template <typename IntegrandFunctor>
void StandardTestKernelTrialIntegral<IntegrandFunctor>::
evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType>& testValues,
        const CollectionOf3dArrays<BasisFunctionType>& trialValues,
        const CollectionOf3dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& quadWeights,
        arma::Mat<ResultType>& result) const
{
    // Evaluate constants

    const int testDofCount = testValues[0].extent(1);
    const int trialDofCount = trialValues[0].extent(1);

    const int pointCount = quadWeights.size();

    // Assert that array dimensions are correct

    for (int i = 0; i < kernelValues.size(); ++i)
        assert(kernelValues[i].extent(2) == pointCount);
    for (int i = 0; i < testValues.size(); ++i)
        assert(testValues[i].extent(2) == pointCount);
    for (int i = 0; i < trialValues.size(); ++i)
        assert(trialValues[i].extent(2) == pointCount);

    // Integrate

    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (int testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (int point = 0; point < pointCount; ++point)
                sum += m_functor.evaluate(
                            testGeomData.const_slice(point),
                            trialGeomData.const_slice(point),
                            testValues.const_slice(testDof, point),
                            trialValues.const_slice(trialDof, point),
                            kernelValues.const_slice(point)) *
                        testGeomData.integrationElements(point) *
                        trialGeomData.integrationElements(point) *
                        quadWeights[point] ;
            result(testDof, trialDof) = sum;
        }
}

} // namespace Fiber
