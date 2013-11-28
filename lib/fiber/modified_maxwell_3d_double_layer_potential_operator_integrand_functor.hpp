// Copyright (C) 2011-2012 by the BEM++ Authors
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

#ifndef fiber_modified_maxwell_3d_double_layer_potential_operator_integrand_functor_hpp
#define fiber_modified_maxwell_3d_double_layer_potential_operator_integrand_functor_hpp

#include "../common/common.hpp"

#include "collection_of_2d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

#include <cassert>
#include <vector>

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class ModifiedMaxwell3dDoubleLayerPotentialOperatorIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(size_t& trialGeomDeps) const {
        // Do nothing
    }

    int resultDimension() const {
        return 3;
    }

    template <template<typename T> class CollectionOf1dSlicesOfConstNdArrays,
              typename TrialValueType>
    void evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& /* trialGeomData */,
            const CollectionOf2dSlicesOfConst4dArrays<KernelType>& kernelValues,
            const CollectionOf1dSlicesOfConstNdArrays<TrialValueType>& trialValues,
            std::vector<ResultType>& result) const {
#       ifndef NDEBUG
        const int dimWorld = 3;
#       endif

        // Assert that there is at least one vector-valued kernel
        assert(kernelValues.size() >= 1);
        assert(kernelValues[0].extent(0) == 3);
        assert(kernelValues[0].extent(1) == 1);

        // Assert that there are is at least one trial transformations
        // (function value)
        assert(trialValues.size() >= 1);
        assert(trialValues[0].extent(0) == 3);

        // Assert that the result vector is three-dimensional
        assert(result.size() == dimWorld);

        result[0] = kernelValues[0](1, 0) * trialValues[0](2) -
                    kernelValues[0](2, 0) * trialValues[0](1);
        result[1] = kernelValues[0](2, 0) * trialValues[0](0) -
                    kernelValues[0](0, 0) * trialValues[0](2);
        result[2] = kernelValues[0](0, 0) * trialValues[0](1) -
                    kernelValues[0](1, 0) * trialValues[0](0);
    }
};

} // namespace Fiber

#endif
