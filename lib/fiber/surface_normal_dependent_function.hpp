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

#ifndef fiber_surface_normal_dependent_function_hpp
#define fiber_surface_normal_dependent_function_hpp

#include "../common/common.hpp"

#include "function.hpp"
#include "geometrical_data.hpp"

namespace Fiber
{

/** \brief Function to be used as a source term, depending on global
  coordinates and on the surface normal.

  The template parameter \p Functor should be a class with the following
  interface:

  class Functor
  {
  public:
      // Type of the function's values (e.g. float or std::complex<double>)
      typedef <implementiation-defined> ValueType;
      typedef ScalarTraits<ValueType>::RealType CoordinateType;

      // Copy constructor
      Functor(const Functor& other);

      // Number of components of the function's arguments ("point" and "normal")
      static const size_t argumentDimension;

      // Number of components of the function's result
      static const size_t resultDimension;

      // Evaluate the function at the point "point", with vector normal to the
      // surface given in the argument "normal", and store result in the array
      // "result".
      // All arrays will be preinitialised to correct dimensions.
      void evaluate(const arma::Col<CoordinateType>& point,
                    const arma::Col<CoordinateType>& normal,
                    arma::Col<ValueType>& result) const;
  };
  */
template <typename Functor>
class SurfaceNormalDependentFunction : public Function<typename Functor::ValueType>
{
public:
    typedef Function<typename Functor::ValueType> Base;
    typedef typename Functor::ValueType ValueType;
    typedef typename Base::CoordinateType CoordinateType;

    SurfaceNormalDependentFunction(const Functor& functor) :
        m_functor(functor) {
    }

    virtual size_t worldDimension() const {
        return m_functor.argumentDimension;
    }

    virtual size_t codomainDimension() const {
        return m_functor.resultDimension;
    }

    virtual void addGeometricalDependencies(size_t& geomDeps) const {
        geomDeps |= GLOBALS | NORMALS;
    }

    virtual void evaluate(const GeometricalData<CoordinateType>& geomData,
                          arma::Mat<ValueType>& result) const {
        const arma::Mat<CoordinateType>& points  = geomData.globals;
        const arma::Mat<CoordinateType>& normals = geomData.normals;

#ifndef NDEBUG
        if (points.n_rows != worldDimension())
            throw std::invalid_argument(
                    "SurfaceNormalDependentFunction::evaluate(): "
                    "incompatible world dimension");
#endif

        const size_t pointCount = points.n_cols;
        result.set_size(codomainDimension(), pointCount);
        for (size_t i = 0; i < pointCount; ++i) {
            arma::Col<ValueType> activeResultColumn = result.unsafe_col(i);
            m_functor.evaluate(points.unsafe_col(i), normals.unsafe_col(i),
                               activeResultColumn);
        }
    }

private:
    Functor m_functor;
};

} // namespace Fiber

#endif
