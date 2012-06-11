#ifndef fiber_simple_scalar_kernel_trial_integrand_functor_hpp
#define fiber_simple_scalar_kernel_trial_integrand_functor_hpp

#include <cassert>
#include "collection_of_2d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class SimpleScalarKernelTrialIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(int& trialGeomDeps) const {
        // do nothing
    }

    int resultDimension() const {
        return 1;
    }

    // It is possible that this function could be generalised to
    // multiple basis transformations or kernels and that the additional
    // loops could be optimised away by the compiler.
    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& /* trialGeomData */,
            const CollectionOf2dSlicesOfConst4dArrays<KernelType>& kernelValues,
            const CollectionOf1dSlicesOfConst2dArrays<ResultType>&
            weightedTransformedTrialValues) const {
        // Assert that there is only a single scalar-valued kernel
        assert(kernelValues.size() == 1);
        assert(kernelValues[0].extent(0) == 1);
        assert(kernelValues[0].extent(1) == 1);

        // Assert that there is only a single weighted trial transformation
        // and that it is scalar
        assert(weightedTransformedTrialValues.size() == 1);
        const int transformationDim = weightedTransformedTrialValues[0].extent(0);
        assert(transformationDim == 1);

        return weightedTransformedTrialValues[0](0) * kernelValues[0](0, 0);
    }
};

} // namespace Fiber

#endif
