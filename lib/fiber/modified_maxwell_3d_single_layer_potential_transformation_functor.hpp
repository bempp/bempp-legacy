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

#ifndef fiber_modified_maxwell_3d_single_layer_potential_transformation_functor_hpp
#define fiber_modified_maxwell_3d_single_layer_potential_transformation_functor_hpp

#include "../common/common.hpp"

#include "basis_transformation_functor_wrappers.hpp"
#include "hdiv_function_value_functor.hpp"
#include "surface_div_3d_functor.hpp"

namespace Fiber
{

/** \brief Functor calculating basis-function transformations necessary for
 *  the implementation of the single-layer potential boundary operator for the
 *  modified Maxwell equations in 3D.
 *
 *  Two following quantities are calculated:
 *  * shape-function value
 *  * shape-function surface div. */
template <typename CoordinateType_>
class ModifiedMaxwell3dSingleLayerPotentialTransformationFunctor :
        public ElementaryBasisTransformationFunctorPairWrapper<
        HdivFunctionValueElementaryFunctor<CoordinateType_>,
        SurfaceDiv3dElementaryFunctor<CoordinateType_> >
{
public:
    typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
