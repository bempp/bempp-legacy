// Copyright (C) 2011-2012 by the Bem++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef fiber_simple_scalar_kernel_trial_integrand_functor_hpp
#define fiber_simple_scalar_kernel_trial_integrand_functor_hpp

#include "../common/common.hpp"

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

    void addGeometricalDependencies(size_t& trialGeomDeps) const {
        // do nothing
    }

    int resultDimension() const {
        return 1;
    }

    // It is possible that this function could be generalised to
    // multiple basis transformations or kernels and that the additional
    // loops could be optimised away by the compiler.
    void evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& /* trialGeomData */,
            const CollectionOf2dSlicesOfConst4dArrays<KernelType>& kernelValues,
            const CollectionOf1dSlicesOfConst2dArrays<ResultType>&
            weightedTransformedTrialValues,
            std::vector<ResultType>& value) const {
        // Assert that there is at least one scalar-valued kernel
        assert(kernelValues.size() >= 1);
        assert(kernelValues[0].extent(0) == 1);
        assert(kernelValues[0].extent(1) == 1);

        // Assert that there is at least one scalar weighted trial
        // transformation
        assert(weightedTransformedTrialValues.size() >= 1);
#ifndef NDEBUG
        const int transformationDim = weightedTransformedTrialValues[0].extent(0);
        assert(transformationDim == 1);
#endif
        assert(value.size() == 1);

        value[0] = weightedTransformedTrialValues[0](0) * kernelValues[0](0, 0);
    }
};

} // namespace Fiber

#endif
