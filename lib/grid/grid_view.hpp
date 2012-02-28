// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_grid_view_hpp
#define bempp_grid_view_hpp

#include "common.hpp"
#include "geometry_type.hpp"
#include "entity_iterator.hpp"

#include <armadillo>
#include <boost/utility/enable_if.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <memory>
#include <stdexcept>

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;
template<int codim> class EntityCache;
class IndexSet;
template<int codim> class Mapper;
class ReverseIndexSet;
class VtkWriter;

/** \brief Abstract wrapper of a grid view */
class GridView
{
public:
    /** Destructor */
    virtual ~GridView() {
    }

    /** \brief The index set */
    virtual const IndexSet& indexSet() const = 0;

    /** \brief The element mapper. */
    virtual const Mapper<0>& elementMapper() const = 0;

    /** \brief Number of entities with codimension \p codim. */
    virtual int entityCount(int codim) const = 0;

    /** \brief Number of entities with geometry type \p type. */
    virtual int entityCount(const GeometryType &type) const = 0;

    /** \brief True if the entity \p e of codimension 0 is contained in this grid view.
       *
       * \note If \p e is not an element of the grid, then
       *       the result of containsEntity() is undefined.
       */
    virtual bool containsEntity(const Entity<0>& e) const = 0;
    /** \brief True if the entity \p e of codimension 1 is contained in this grid view.

      \overload
    */
    virtual bool containsEntity(const Entity<1>& e) const = 0;
    /** \brief True if the entity \p e of codimension 2 is contained in this grid view.

      \overload
    */
    virtual bool containsEntity(const Entity<2>& e) const = 0;
    /** \brief True if the entity \p e of codimension 3 is contained in this grid view.

      \overload
    */
    virtual bool containsEntity(const Entity<3>& e) const = 0;

    /** \brief Iterator over entities of codimension \p codim contained in this view. */
    // Default implementation; specialisations for potentially allowed codimensions follow
    // after class declaration.
    template<int codim>
    std::auto_ptr<EntityIterator<codim> > entityIterator() const {
        throw std::logic_error("GridView::entityIterator(): invalid entity codimension");
    }

    /** \brief Get raw data describing the geometry of all codim-0 entities
      contained in this grid view.

      This method is mainly intended for use in the OpenCL implementation.

      \param[out] vertices
        On output, a 2D array whose (i,j)th element is the ith
        coordinate of the vertex of index j.

      \param[out] elementCorners
        On output, a 2D array whose (i,j)th element is the index of the ith
        corner of jth codim-0 entity, or -1 if this entity has less than
        i-1 corners.

      \note For isoparametric elements we will likely need to add a third
        parameter to contain "arbitrary" auxiliary data.
      */
    virtual void getRawElementData(arma::Mat<ctype>& vertices,
                                   arma::Mat<int>& elementCorners,
                                   arma::Mat<char>& auxData) const = 0;

    /** \brief The reverse index set.

      Note that this object is *not* updated when the grid is adapted. In that
      case you need to create a new grid view and obtain a new reference to the
      reverse index set.

      \internal Actual creation of the reverse index set is delayed to the
      first call to this method, since the reverse index set needs to store the
      entity pointers to all entities, and hence takes up a lot of memory.
    */
    virtual const ReverseIndexSet& reverseIndexSet() const = 0;

    /** \brief Create a VtkWriter for this grid view.

      \param dm Data mode (conforming or nonconforming; see the documentation of Dune::VTK::DataMode for details). */
    virtual std::auto_ptr<VtkWriter> vtkWriter(Dune::VTK::DataMode dm=Dune::VTK::conforming) const = 0;

    // Deferred for later implementation:
    // * Iteration over neighbours: Dune methods ibegin() and iend().

private:
    /** \brief Iterator over entities of codimension 0 contained in this view. */
    virtual std::auto_ptr<EntityIterator<0> > entityCodim0Iterator() const = 0;
    /** \brief Iterator over entities of codimension 1 contained in this view. */
    virtual std::auto_ptr<EntityIterator<1> > entityCodim1Iterator() const = 0;
    /** \brief Iterator over entities of codimension 2 contained in this view. */
    virtual std::auto_ptr<EntityIterator<2> > entityCodim2Iterator() const = 0;
    /** \brief Iterator over entities of codimension 3 contained in this view. */
    virtual std::auto_ptr<EntityIterator<3> > entityCodim3Iterator() const = 0;
};

template<>
inline std::auto_ptr<EntityIterator<0> > GridView::entityIterator<0>() const
{
    return entityCodim0Iterator();
}
template<>
inline std::auto_ptr<EntityIterator<1> > GridView::entityIterator<1>() const
{
    return entityCodim1Iterator();
}
template<>
inline std::auto_ptr<EntityIterator<2> > GridView::entityIterator<2>() const
{
    return entityCodim2Iterator();
}
template<>
inline std::auto_ptr<EntityIterator<3> > GridView::entityIterator<3>() const
{
    return entityCodim3Iterator() ;
}

} // namespace Bempp

#endif
