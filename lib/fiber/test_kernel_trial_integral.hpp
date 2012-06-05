#ifndef fiber_test_kernel_trial_integral_hpp
#define fiber_test_kernel_trial_integral_hpp

namespace Fiber
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral
{
public:
    virtual void evaluateWithTensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            const Array3dCollection<BasisFunctionType>& testExpressions,
            const Array3dCollection<BasisFunctionType>& trialExpressions,
            const Array4dCollection<KernelType>& kernels,
            const std::vector<CoordinateType>& testQuadWeights,
            const std::vector<CoordinateType>& trialQuadWeights,
            const Slice2d<ResultType>& result) const = 0;

    virtual void evaluateWithNontensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            const Array3dCollection<BasisFunctionType>& testExpressions,
            const Array3dCollection<BasisFunctionType>& trialExpressions,
            const Array3dCollection<KernelType>& kernels,
            const std::vector<CoordinateType>& quadWeights,
            const Slice2d<ResultType>& result) const = 0;
};

} // namespace Fiber

#endif
