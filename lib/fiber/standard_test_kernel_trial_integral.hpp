#ifndef fiber_standard_test_kernel_trial_integral_hpp
#define fiber_standard_test_kernel_trial_integral_hpp

#include "test_kernel_trial_integral.hpp"

namespace Fiber
{

template <typename IntegrandFunctor>
class StandardTestKernelTrialIntegral :
        public TestKernelTrialIntegral<
        typename IntegrandFunctor::BasisFunctionType,
        typename IntegrandFunctor::KernelType,
        typename IntegrandFunctor::ResultType>
{
public:
    virtual void evaluateWithTensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const Array3dCollection<BasisFunctionType>& testExpressions,
            const Array4dCollection<KernelType>& kernels,
            const Array3dCollection<BasisFunctionType>& trialExpressions,
            const GeometricalData<CoordinateType>& trialGeomData,
            const std::vector<CoordinateType>& quadWeights,
            arma::Cube<ResultType>& result) const;

    virtual void evaluateWithNontensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const Array3dCollection<BasisFunctionType>& testExpressions,
            const Array3dCollection<KernelType>& kernels,
            const Array3dCollection<BasisFunctionType>& trialExpressions,
            const GeometricalData<CoordinateType>& trialGeomData,
            const std::vector<CoordinateType>& quadWeights,
            arma::Cube<ResultType>& result) const;
};

} // namespace Fiber

#endif
