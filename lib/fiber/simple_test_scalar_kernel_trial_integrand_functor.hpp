#ifndef fiber_simple_test_scalar_kernel_trial_integrand_functor_hpp
#define fiber_simple_test_scalar_kernel_trial_integrand_functor_hpp

#include <cassert>
#include "collection_of_3d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class SimpleTestScalarKernelTrialIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(int& testGeomDeps, int& trialGeomDeps) const {
        // do nothing
    }

    // It is possible that this function could be generalised to
    // multiple basis transformations or kernels and that the additional
    // loops could be optimised away by the compiler.
    template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& /* testGeomData */,
            const ConstGeometricalDataSlice<CoordinateType>& /* trialGeomData */,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& testValues,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& trialValues,
            const CollectionOf2dSlicesOfConstNdArrays<KernelType>& kernelValues) const {
        // Assert that there is only a single scalar-valued kernel
        assert(kernelValues.size() == 1);
        assert(kernelValues[0].extent(0) == 1);
        assert(kernelValues[0].extent(1) == 1);

        // Assert that there is only a single test and trial transformation
        // and that their dimensions agree
        assert(testValues.size() == 1);
        assert(trialValues.size() == 1);
        assert(testValues[0].extent(0) == trialValues[0].extent(0));

        const int transformationDim = testValues[0].extent(0);

        ResultType result = 0.;
        for (int dim = 0; dim < transformationDim; ++dim)
            result += conjugate(testValues[0](dim)) *
                    trialValues[0](dim);
        result *= kernelValues[0](0, 0);
        return result;
    }
};

} // namespace Fiber

#endif
