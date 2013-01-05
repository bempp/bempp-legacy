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

#ifndef fiber_modified_maxwell_3d_single_layer_potential_integrand_functor_hpp
#define fiber_modified_maxwell_3d_single_layer_potential_integrand_functor_hpp

#include "../common/common.hpp"

#include <cassert>
#include "collection_of_3d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"

namespace Fiber
{

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class ModifiedMaxwell3dSingleLayerPotentialIntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const {
        testGeomDeps |= NORMALS;
        trialGeomDeps |= NORMALS;
    }

    template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& testTransfValues,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>& trialTransfValues,
            const CollectionOf2dSlicesOfConstNdArrays<KernelType>& kernelValues) const {
        const int dimWorld = 3;

        // Assert that there are at least two scalar-valued kernels
        assert(kernelValues.size() >= 2);
        assert(kernelValues[0].extent(0) == 1);
        assert(kernelValues[0].extent(1) == 1);
        assert(kernelValues[1].extent(0) == 1);
        assert(kernelValues[1].extent(1) == 1);

        // Assert that there are at least two test and trial transformations
        // (function value and surface div) of correct dimensions
        assert(testTransfValues.size() >= 2);
        assert(trialTransfValues.size() >= 2);
        _1dSliceOfConst3dArray<BasisFunctionType> testValues = testTransfValues[0];
        _1dSliceOfConst3dArray<BasisFunctionType> trialValues = trialTransfValues[0];
        _1dSliceOfConst3dArray<BasisFunctionType> testSurfaceDivs = testTransfValues[1];
        _1dSliceOfConst3dArray<BasisFunctionType> trialSurfaceDivs = trialTransfValues[1];
        assert(testValues.extent(0) == 3);
        assert(trialValues.extent(0) == 3);
        assert(testSurfaceDivs.extent(0) == 1);
        assert(trialSurfaceDivs.extent(0) == 1);

        // Let K_0(x, y) = K(x, y) and K_1(x, y) = K(x, y) / kappa^2.
        // Return
        // K_0(x, y) u*(x) . v(y) - K_1(x, y) div u*(x) div v(y)

        ResultType term_0 = 0.;
        for (int dim = 0; dim < dimWorld; ++dim)
            term_0 += conjugate(testValues(dim)) * trialValues(dim);
        term_0 *= kernelValues[0](0, 0);
        ResultType term_1 = conjugate(testSurfaceDivs(0)) * trialSurfaceDivs(0) *
            kernelValues[1](0, 0);
        return term_0 - term_1;
    }
};

} // namespace Fiber

#endif
