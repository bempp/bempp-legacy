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

#ifndef bempp_surface_normal_independent_function_hpp
#define bempp_surface_normal_independent_function_hpp

#include "../common/common.hpp"

#include "../fiber/surface_normal_independent_function.hpp"

namespace Bempp {

using Fiber::SurfaceNormalIndependentFunction;

/** \ingroup assembly_functions
  \brief Use a functor depending on global coordinates only (independent from
  the surface normal or the element index) to construct an instance of
  the Function class.

  The template parameter \p Functor should be a class with the following
  interface:

  \code
  class Functor
  {
  public:
      // Type of the function's values (float, double, std::complex<float>
      // or std::complex<double>)
      typedef <implementiation-defined> ValueType;
      typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

      // Number of components of the function's argument
      int argumentDimension() const;

      // Number of components of the function's result
      int resultDimension() const;

      // Evaluate the function at the point "point" and store result in
      // the array "result". The "point" and "result" arrays will be
      // preinitialized to correct dimensions.
      void evaluate(const Vector<CoordinateType>& point,
                    Vector<ValueType>& result) const;
  };
  \endcode

  The constructed Function object can subsequently be passed into a constructor
  of the GridFunction class. */
template <typename Functor>
inline SurfaceNormalIndependentFunction<Functor>
surfaceNormalIndependentFunction(const Functor &functor) {
  return SurfaceNormalIndependentFunction<Functor>(functor);
}

} // namespace Bempp

#endif
