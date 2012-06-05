#ifndef fiber_integrand_hpp
#define fiber_integrand_hpp

namespace Fiber
{

template <typename ResultType>
class TestKernelTrialIntegral
{
public:
    virtual void evaluateWithTensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const Array3dCollection<BasisFunctionType>& testExpressions,
            const Array4dCollection<KernelType>& kernels,
            const Array3dCollection<BasisFunctionType>& trialExpressions,
            const GeometricalData<CoordinateType>& trialGeomData,
            const std::vector<CoordinateType>& quadWeights,
            arma::Cube<ResultType>& result) const = 0;

    virtual void evaluateWithNontensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const Array3dCollection<BasisFunctionType>& testExpressions,
            const Array3dCollection<KernelType>& kernels,
            const Array3dCollection<BasisFunctionType>& trialExpressions,
            const GeometricalData<CoordinateType>& trialGeomData,
            const std::vector<CoordinateType>& quadWeights,
            arma::Cube<ResultType>& result) const = 0;
};

} // namespace Fiber

#endif
