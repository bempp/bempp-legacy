// Copyright (C) 2011-2012 by the Fiber Authors
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

#ifndef fiber_ordinary_function_hpp
#define fiber_ordinary_function_hpp

#include "function.hpp"
#include "geometrical_data.hpp"

namespace Fiber
{

/** \brief Function defined via a user-supplied functor, depending only on
    global coordinates.

  The template parameter \p Functor should be a class with the following
  interface:

  class Functor
  {
  public:
      // Type of the function's values (float, double, std::complex<float>
      // or std::complex<double>)
      typedef <implementiation-defined> ValueType;
      typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

      // Copy constructor
      Functor(const Functor& other);

      // Number of components of the function's argument
      static const int argumentDimension;

      // Number of components of the function's result
      static const int resultDimension;

      // Evaluate the function at the point "point" and store result in
      // the array "result"
      // The "point" and "result" arrays will be preinitialised to correct
      // dimensions.
      void evaluate(const arma::Col<CoordinateType>& point,
                    arma::Col<ValueType>& result) const;
  };
  */
template <typename Functor>
class OrdinaryFunction : public Function<typename Functor::ValueType>
{
public:
    typedef Function<typename Functor::ValueType> Base;
    typedef typename Functor::ValueType ValueType;
    typedef typename Base::CoordinateType CoordinateType;

    OrdinaryFunction(const Functor& functor) :
        m_functor(functor) {
    }

    virtual int worldDimension() const {
        return m_functor.argumentDimension;
    }

    virtual int codomainDimension() const {
        return m_functor.resultDimension;
    }

    virtual void addGeometricalDependencies(int& geomDeps) const {
        geomDeps |= GLOBALS;
    }

    virtual void evaluate(const GeometricalData<CoordinateType>& geomData,
                          arma::Mat<ValueType>& result) const {
        const arma::Mat<CoordinateType>& points = geomData.globals;

#ifndef NDEBUG
        if (points.n_rows != worldDimension())
            throw std::invalid_argument("OrdinaryScalarFunction::evaluate(): "
                                        "incompatible world dimension");
#endif

        const int pointCount = points.n_cols;
        result.set_size(codomainDimension(), pointCount);
        for (int i = 0; i < pointCount; ++i)
        {
            arma::Col<ValueType> activeResultColumn = result.unsafe_col(i);
            m_functor.evaluate(points.unsafe_col(i), activeResultColumn);
        }
    }

private:
    Functor m_functor;
};

} // namespace Fiber

#endif
