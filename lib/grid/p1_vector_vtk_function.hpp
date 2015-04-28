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

#ifndef bempp_p1_vector_vtk_function_hpp
#define bempp_p1_vector_vtk_function_hpp

#include "../common/common.hpp"

#include <dune/grid/io/file/vtk/function.hh>

namespace Bempp {

/**
 * \ingroup grid_internal
 * \brief Take a vector and interpret it as point data for the VTKWriter
 *
 * This class turns a generic vector containing point data into a
 * VTKFunction.  The vector must allow read access to the data via
 * operator[]() and store the data in the order given by
 * MultipleCodimMultipleGeomTypeMapper with a layout class that allows only
 * vertices.  Also, it must support the method size().
 *
 * As opposed to a <tt>Dune::P1VTKFunction</tt>, this class represents all
 * components of a vector field at once.
 *
 * \note While this function may be evaluated anywhere on a given grid
 *       element, it does not interpolate between the corners of the element
 *       -- instead, it returns the value at the nearest corner (as
 *       determined in local coordinates).
 *
 * \tparam GV Type of GridView the vector applies to.
 * \tparam V  Type of vector.
 */
template <typename GV, typename V>
class P1VectorVTKFunction : public Dune::VTKFunction<GV> {
  //! Base class
  typedef Dune::VTKFunction<GV> Base;
  //! Mapper for vertices
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV, Dune::MCMGVertexLayout>
  Mapper;

  //! store a reference to the vector
  const V &v;
  //! name of this function
  std::string s;
  //! number of components of the field stored in the vector
  int ncomps_;
  //! mapper used to map elements to indices
  Mapper mapper;

public:
  typedef typename Base::Entity Entity;
  typedef typename Base::ctype ctype;
  using Base::dim;

  //! return number of components
  virtual int ncomps() const { return ncomps_; }

  //! evaluate
  virtual double evaluate(int comp, const Entity &e,
                          const Dune::FieldVector<ctype, dim> &xi) const {
    double min = 1E100;
    int imin = -1;
    Dune::GeometryType gt = e.type();
    for (int i = 0; i < e.template count<dim>(); ++i) {
      Dune::FieldVector<ctype, dim> local =
          Dune::ReferenceElements<ctype, dim>::general(gt)
              .position(i, dim);
      local -= xi;
      if (local.infinity_norm() < min) {
        min = local.infinity_norm();
        imin = i;
      }
    }
    return v(mapper.map(e, imin, dim) * ncomps_ + comp);
  }

  //! get name
  virtual std::string name() const { return s; }

  /**
   * \brief Construct from a vector and a name.
   * \param gv     GridView to operate on (used to instantiate a
   *               MultipleCodimMultipleGeomTypeMapper, otherwise no
   *               reference or copy is stored).  Note that this must be the
   *               GridView the vector applies to as well as the GridView
   *               later used by the VTKWriter -- i.e. we do not implicitly
   *               restrict or prolongate the data.
   * \param v_     Reference to the vector holding the data.  The reference
   *               is stored internally and must be valid for as long as
   *               this functions evaluate method is used.
   * \param s_     Name of this function in the VTK file.
   * \param ncomps Number of components of the field represented by the
   *               vector. */
  P1VectorVTKFunction(const GV &gv, const V &v_, const std::string &s_,
                      int ncomps = 1)
      : v(v_), s(s_), ncomps_(ncomps), mapper(gv) {
    if (v.size() != (unsigned int)(mapper.size() * ncomps_))
      DUNE_THROW(Dune::IOError, "P1VectorVTKFunction: size mismatch");
  }
};

} // namespace Bempp

#endif
