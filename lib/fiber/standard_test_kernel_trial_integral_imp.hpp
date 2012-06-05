#include "standard_test_kernel_trial_integral.hpp"

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class IntegrandFunctor
{
public:
    typedef ResultType_ ResultType;
    ResultType evaluate(
            const GeometricalDataSlice<CoordinateType>& testGeomData,
            const GeometricalDataSlice<CoordinateType>& trialGeomData,
            const Slice1dCollection<BasisFunctionType>& testExpressions,
            const Slice1dCollection<BasisFunctionType>& trialExpressions,
            const Slice2dCollection<KernelType>& kernels);
};

template <typename IntegrandFunctor>
void StandardTestKernelTrialIntegral<IntegrandFunctor>::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        const Array3dCollection<BasisFunctionType>& testExpressions,
        const Array3dCollection<BasisFunctionType>& trialExpressions,
        const Array4dCollection<KernelType>& kernels,
        const std::vector<CoordinateType>& testQuadWeights,
        const std::vector<CoordinateType>& trialQuadWeights,
        const Slice2d<ResultType>& result) const
{
    // TODO: In debug mode, check that all the array dimensions are consistent
    const int testDofCount = testExpressions.extent(1);
    const int trialDofCount = trialExpressions.extent(1);

    const int testPointCount = testQuadWeights.size();
    const int trialPointCount = trialQuadWeights.size();

    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (int testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
                for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
                    sum += m_integrand.evaluate(
                                testGeomData.slice(testPoint),
                                testExpressions.slice(testDof, testPoint),
                                kernels.slice(testPoint, trialPoint),
                                trialExpressions.slice(trialDof, trialPoint),
                                trialGeomData.slice(trialPoint)) *
                            testGeomData.integrationElements(testPoint) *
                            trialGeomData.integrationElements(trialPoint) *
                            testQuadWeights[testPoint] *
                            trialQuadWeights[trialPoint];
            result(testDof, trialDof) = sum;
        }
}

template <typename IntegrandFunctor>
void StandardTestKernelTrialIntegral<IntegrandFunctor>::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        const Array3dCollection<BasisFunctionType>& testExpressions,
        const Array3dCollection<BasisFunctionType>& trialExpressions,
        const Array4dCollection<KernelType>& kernels,
        const std::vector<CoordinateType>& quadWeights,
        const Slice2d<ResultType>& result) const
{
    // TODO: In debug mode, check that all the array dimensions are consistent
    const int testDofCount = testExpressions.extent(1);
    const int trialDofCount = trialExpressions.extent(1);

    const int pointCount = quadWeights.size();

    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (int testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (int point = 0; point < pointCount; ++point)
                sum += m_integrand.evaluate(
                            testGeomData.slice(point),
                            testExpressions.slice(testDof, point),
                            kernels.slice(point),
                            trialExpressions.slice(trialDof, point),
                            trialGeomData.slice(point)) *
                        testGeomData.integrationElements(point) *
                        trialGeomData.integrationElements(point) *
                        quadWeights[point] ;
            result(testDof, trialDof) = sum;
        }
}

} // namespace Fiber

#endif
