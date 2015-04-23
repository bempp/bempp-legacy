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

#ifndef fiber_surface_normal_and_domain_index_dependent_function_hpp
#define fiber_surface_normal_and_domain_index_dependent_function_hpp

#include "../common/common.hpp"

#include "function.hpp"
#include "geometrical_data.hpp"
#include "types.hpp"

namespace Fiber {

/** \brief %Function intended to be evaluated on a boundary-element grid,
  defined via a user-supplied functor depending on the global coordinates
  of a point lying on the grid, the local unit vector normal to the grid and the
  index of the domain containing the point.

  The template parameter \p Functor should be a class implementing the following
  interface:

  \code
  class Functor
  {
  public:
      // Type of the function's values (e.g. float or std::complex<double>)
      typedef <implementiation-defined> ValueType;
      typedef ScalarTraits<ValueType>::RealType CoordinateType;

      // Number of components of the function's arguments ("point" and "normal")
      int argumentDimension() const;

      // Number of components of the function's result
      int resultDimension() const;

      // Evaluate the function at the point "point" lying on an element from
      // domain "domainIndex", with "normal" the local unit normal vector, and
      // store result in the array "result". All arrays will be preinitialised
      // to correct dimensions.
      void evaluate(const Eigen::Ref<Vector<CoordinateType>>& point,
                    const Eigen::Ref<Vector<CoordinateType>>& normal,
                    int domainIndex,
                    Eigen::Ref<Vector<ValueType>> result) const;
  };
  \endcode
*/
template <typename Functor>
class SurfaceNormalAndDomainIndexDependentFunction
    : public Function<typename Functor::ValueType> {
public:
  typedef Function<typename Functor::ValueType> Base;
  typedef typename Base::ValueType ValueType;
  typedef typename Base::CoordinateType CoordinateType;

  SurfaceNormalAndDomainIndexDependentFunction(const Functor &functor)
      : m_functor(functor) {}

  virtual int worldDimension() const { return m_functor.argumentDimension(); }

  virtual int codomainDimension() const { return m_functor.resultDimension(); }

  virtual void addGeometricalDependencies(size_t &geomDeps) const {
    geomDeps |= GLOBALS | NORMALS | DOMAIN_INDEX;
  }

  virtual void evaluate(const GeometricalData<CoordinateType> &geomData,
                        Matrix<ValueType> &result) const {
    const Matrix<CoordinateType> &points = geomData.globals;
    const Matrix<CoordinateType> &normals = geomData.normals;

#ifndef NDEBUG
    if ((int)points.rows() != worldDimension() ||
        (int)points.rows() != worldDimension())
      throw std::invalid_argument(
          "SurfaceNormalAndDomainIndexDependentFunction::evaluate(): "
          "incompatible world dimension");
#endif

    const size_t pointCount = points.cols();
    result.resize(codomainDimension(), pointCount);
    for (size_t i = 0; i < pointCount; ++i) {
      m_functor.evaluate(
              Eigen::Ref<Vector<CoordinateType>>(const_cast<Matrix<CoordinateType>&>(points).col(i)), 
              Eigen::Ref<Vector<CoordinateType>>(const_cast<Matrix<CoordinateType>&>(normals).col(i)),
              geomData.domainIndex, 
              Eigen::Ref<Vector<ValueType>>(result.col(i)));
    }
  }

private:
  Functor m_functor;
};

} // namespace Fiber

#endif
