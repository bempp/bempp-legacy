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

#ifndef bempp_p0_vector_vtk_function_hpp
#define bempp_p0_vector_vtk_function_hpp

#include "../common/common.hpp"

#include <dune/grid/io/file/vtk/function.hh>

namespace Bempp {

/**
 * \ingroup grid_internal
 * \brief Take a vector and interpret it as cell data for the VTKWriter.
 *
 * This class turns a generic vector containing cell data into a
 * VTKFunction.  The vector must allow read access to the data via
 * operator[]() and store the data in the order given by
 * MultipleCodimMultipleGeomTypeMapper with a layout class that allows only
 * elements.  Also, it must support the method size().
 *
 * As opposed to a <tt>Dune::P0VTKFunction</tt>, this class represents all
 * components of a vector field at once.
 *
 * \tparam GV Type of GridView the vector applies to.
 * \tparam V  Type of vector.
 */
template <typename GV, typename V>
class P0VectorVTKFunction : public Dune::VTKFunction<GV> {
  //! Base class
  typedef Dune::VTKFunction<GV> Base;
  //! Mapper for elements
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV, Dune::MCMGElementLayout>
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
    return v(mapper.map(e) * ncomps_ + comp);
  }

  //! get name
  virtual std::string name() const { return s; }

  /**
   * \brief Construct from a vector and a name.
   *
   * \param gv     GridView to operate on (used to instantiate a
   *               MultipleCodimMultipleGeomeTypeMapper, otherwise no
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
  P0VectorVTKFunction(const GV &gv, const V &v_, const std::string &s_,
                      int ncomps = 1)
      : v(v_), s(s_), ncomps_(ncomps), mapper(gv) {
    if (v.size() != (unsigned int)(mapper.size() * ncomps_))
      DUNE_THROW(Dune::IOError, "VectorP0VTKFunction: size mismatch");
  }
};

} // namespace Bempp

#endif
