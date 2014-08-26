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

#ifndef fiber_modified_helmholtz_3d_hypersingular_transformation_functor_2_hpp
#define fiber_modified_helmholtz_3d_hypersingular_transformation_functor_2_hpp

#include "../common/common.hpp"

#include "shape_transformation_functor_wrappers.hpp"
#include "scalar_function_value_times_normal_functor.hpp"
#include "surface_curl_3d_functor.hpp"

namespace Fiber {

/** \brief Functor calculating shape-function transformations necessary for
 *  the implementation of the hypersingular operator for the modified Helmholtz
 *  equation in 3D.
 *
 *  Two following quantities are calculated:
 *  * basis-function surface curl
 *  * basis-function value times surface normal. */
template <typename CoordinateType_>
class ModifiedHelmholtz3dHypersingularTransformationFunctor2
    : public ElementaryShapeTransformationFunctorPairWrapper<
          SurfaceCurl3dElementaryFunctor<CoordinateType_>,
          ScalarFunctionValueTimesNormalElementaryFunctor<CoordinateType_>> {
public:
  typedef CoordinateType_ CoordinateType;
};

} // namespace Fiber

#endif
