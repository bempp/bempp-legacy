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
#ifndef bempp_grid_view_hpp
#define bempp_grid_view_hpp

#include "../common/common.hpp"

#include "dune.hpp"
#include "entity_iterator.hpp"

#include "../common/armadillo_fwd.hpp"
#include <boost/utility/enable_if.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <memory>
#include <stdexcept>

namespace Bempp {

/** \cond FORWARD_DECL */
template <int codim> class Entity;
template <int codim> class EntityCache;
class IndexSet;
class Mapper;
class ReverseElementMapper;
class VtkWriter;
/** \endcond */

/** \ingroup grid
    \brief Abstract wrapper of a grid view */
class GridView {
public:
  /** Destructor */
  virtual ~GridView() {}

  /** \brief Dimension of the grid. */
  virtual int dim() const = 0;

  /** \brief Dimension of the space containing the grid. */
  virtual int dimWorld() const = 0;

  /** \brief The index set */
  virtual const IndexSet &indexSet() const = 0;

  /** \brief The element mapper. */
  virtual const Mapper &elementMapper() const = 0;

  /** \brief Number of entities with codimension \p codim. */
  virtual size_t entityCount(int codim) const = 0;

  /** \brief Number of entities with geometry type \p type. */
  virtual size_t entityCount(const GeometryType &type) const = 0;

  /** \brief True if the entity \p e of codimension 0 is contained in this grid
    *view.
     *
     * \note If \p e is not an element of the grid, then
     *       the result of containsEntity() is undefined.
     */
  virtual bool containsEntity(const Entity<0> &e) const = 0;
  /** \brief True if the entity \p e of codimension 1 is contained in this grid
    view.

    \overload
  */
  virtual bool containsEntity(const Entity<1> &e) const = 0;
  /** \brief True if the entity \p e of codimension 2 is contained in this grid
    view.

    \overload
  */
  virtual bool containsEntity(const Entity<2> &e) const = 0;
  /** \brief True if the entity \p e of codimension 3 is contained in this grid
    view.

    \overload
  */
  virtual bool containsEntity(const Entity<3> &e) const = 0;

  /** \brief Iterator over entities of codimension \p codim contained in this
   * view. */
  // Default implementation; specialisations for potentially allowed
  // codimensions follow
  // after class declaration.
  template <int codim>
  std::unique_ptr<EntityIterator<codim>> entityIterator() const {
    throw std::logic_error(
        "GridView::entityIterator(): invalid entity codimension");
  }

  /** \brief Get raw data describing the geometry of all codim-0 entities
    contained in this grid view.

    \param[out] vertices
      On output, a 2D array whose (i,j)th element is the ith
      coordinate of the vertex of index j.

    \param[out] elementCorners
      On output, a 2D array whose (i,j)th element is the index of the ith
      corner of jth codim-0 entity, or -1 if this entity has less than
      i-1 corners.

    \param[out] auxData
      Auxiliary data, to be used in future to describe the geometry of
      curvilinear elements. Currently this argument is always set to an empty
      matrix.
  */

  /** \overload */
  void getRawElementData(arma::Mat<float> &vertices,
                         arma::Mat<int> &elementCorners,
                         arma::Mat<char> &auxData) const;
  /** \overload */
  void getRawElementData(arma::Mat<double> &vertices,
                         arma::Mat<int> &elementCorners,
                         arma::Mat<char> &auxData) const;

  /** \brief Get raw data describing the geometry of all codim-0 entities
    contained in this grid view.

    This overload takes an additional argument, \p domainIndices, which on
    output contains indices of the grid domains to which the corresponding
    elements from \p elementCorners belong. */
  void getRawElementData(arma::Mat<double> &vertices,
                         arma::Mat<int> &elementCorners,
                         arma::Mat<char> &auxData,
                         std::vector<int> &domainIndices) const;
  /** \overload */
  void getRawElementData(arma::Mat<float> &vertices,
                         arma::Mat<int> &elementCorners,
                         arma::Mat<char> &auxData,
                         std::vector<int> &domainIndices) const;

  /** \brief Mapping from codim-0 entity index to entity pointer.

    Note that this object is *not* updated when the grid is adapted. In that
    case you need to create a new grid view and obtain a new reference to the
    reverse element mapper.

    \internal Actual creation of the reverse element mapper is delayed to the
    first call to this method, since the object needs to store the
    entity pointers to all entities, and hence may take up a lot of memory.
  */
  virtual const ReverseElementMapper &reverseElementMapper() const = 0;

  /** \brief Create a VtkWriter for this grid view.

    \param dm Data mode (conforming or nonconforming; see the documentation of
    Dune::VTK::DataMode for details). */
  virtual std::unique_ptr<VtkWriter>
  vtkWriter(Dune::VTK::DataMode dm = Dune::VTK::conforming) const = 0;

  // Deferred for later implementation:
  // * Iteration over neighbours: Dune methods ibegin() and iend().

private:
  virtual void getRawElementDataDoubleImpl(
      arma::Mat<double> &vertices, arma::Mat<int> &elementCorners,
      arma::Mat<char> &auxData, std::vector<int> *domainIndices) const = 0;
  virtual void getRawElementDataFloatImpl(
      arma::Mat<float> &vertices, arma::Mat<int> &elementCorners,
      arma::Mat<char> &auxData, std::vector<int> *domainIndices) const = 0;

  /** \brief Iterator over entities of codimension 0 contained in this view. */
  virtual std::unique_ptr<EntityIterator<0>> entityCodim0Iterator() const = 0;
  /** \brief Iterator over entities of codimension 1 contained in this view. */
  virtual std::unique_ptr<EntityIterator<1>> entityCodim1Iterator() const = 0;
  /** \brief Iterator over entities of codimension 2 contained in this view. */
  virtual std::unique_ptr<EntityIterator<2>> entityCodim2Iterator() const = 0;
  /** \brief Iterator over entities of codimension 3 contained in this view. */
  virtual std::unique_ptr<EntityIterator<3>> entityCodim3Iterator() const = 0;
};

inline void GridView::getRawElementData(arma::Mat<double> &vertices,
                                        arma::Mat<int> &elementCorners,
                                        arma::Mat<char> &auxData) const {
  getRawElementDataDoubleImpl(vertices, elementCorners, auxData, 0);
}

inline void GridView::getRawElementData(arma::Mat<float> &vertices,
                                        arma::Mat<int> &elementCorners,
                                        arma::Mat<char> &auxData) const {
  getRawElementDataFloatImpl(vertices, elementCorners, auxData, 0);
}

inline void GridView::getRawElementData(arma::Mat<double> &vertices,
                                        arma::Mat<int> &elementCorners,
                                        arma::Mat<char> &auxData,
                                        std::vector<int> &domainIndices) const {
  getRawElementDataDoubleImpl(vertices, elementCorners, auxData,
                              &domainIndices);
}

inline void GridView::getRawElementData(arma::Mat<float> &vertices,
                                        arma::Mat<int> &elementCorners,
                                        arma::Mat<char> &auxData,
                                        std::vector<int> &domainIndices) const {
  getRawElementDataFloatImpl(vertices, elementCorners, auxData, &domainIndices);
}

template <>
inline std::unique_ptr<EntityIterator<0>> GridView::entityIterator<0>() const {
  return entityCodim0Iterator();
}
template <>
inline std::unique_ptr<EntityIterator<1>> GridView::entityIterator<1>() const {
  return entityCodim1Iterator();
}
template <>
inline std::unique_ptr<EntityIterator<2>> GridView::entityIterator<2>() const {
  return entityCodim2Iterator();
}
template <>
inline std::unique_ptr<EntityIterator<3>> GridView::entityIterator<3>() const {
  return entityCodim3Iterator();
}

} // namespace Bempp

#endif
