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

#ifndef fiber_surface_grad_3d_functor_hpp
#define fiber_surface_grad_3d_functor_hpp

#include "../common/common.hpp"

#include "basis_data.hpp"
#include "geometrical_data.hpp"
#include "collection_of_3d_arrays.hpp"
#include "shape_transformation_functor_wrappers.hpp"

#include <boost/array.hpp>

namespace Fiber
{

template <typename CoordinateType_>
class SurfaceGrad3dElementaryFunctor
{
public:
    typedef CoordinateType_ CoordinateType;

    int argumentDimension() const { return 1; }
    int resultDimension() const { return 3; }

    void addDependencies(size_t& basisDeps, size_t& geomDeps) const {
        basisDeps |= DERIVATIVES;
        geomDeps |= JACOBIAN_INVERSES_TRANSPOSED;
    }

    template <typename ValueType>
    void evaluate(
            const ConstBasisDataSlice<ValueType>& basisData,
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            _1dSliceOf3dArray<ValueType>& result) const {
        assert(basisData.componentCount() == 1);
        assert(geomData.dimWorld() == 3);

        // result := gradient of the shape function extended outside
        // the surface so that its normal derivative on the surf. is zero
        result(0) = basisData.derivatives(0, 0) * geomData.jacobianInverseTransposed(0, 0) +
                basisData.derivatives(0, 1) * geomData.jacobianInverseTransposed(0, 1);
        result(1) = basisData.derivatives(0, 0) * geomData.jacobianInverseTransposed(1, 0) +
                basisData.derivatives(0, 1) * geomData.jacobianInverseTransposed(1, 1);
        result(2) = basisData.derivatives(0, 0) * geomData.jacobianInverseTransposed(2, 0) +
                basisData.derivatives(0, 1) * geomData.jacobianInverseTransposed(2, 1);
    }
};

// Note: in C++11 we'll be able to make a "template typedef", or more precisely
// a using declaration, instead of this spurious inheritance
/** \ingroup functors
 *  \brief Functor calculating the surface gradient of a scalar field in 3D. */
template <typename CoordinateType_>
class SurfaceGrad3dFunctor :
        public ElementaryShapeTransformationFunctorWrapper<
        SurfaceGrad3dElementaryFunctor<CoordinateType_> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
