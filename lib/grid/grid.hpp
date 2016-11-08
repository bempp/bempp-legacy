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

#ifndef bempp_grid_hpp
#define bempp_grid_hpp

/** \file . */

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"
#include "grid_parameters.hpp"
#include "entity.hpp"

#include <cstddef> // size_t
#include <memory>
#include <vector>
#include <tbb/mutex.h>
#include <boost/signals2/signal.hpp>
#include <functional>

namespace Bempp {

/** \cond FORWARD_DECL */
template <int codim> class Entity;
template <typename BasisFunctionType> class Space;
class GeometryFactory;
class GridView;
class IdSet;
template <typename CoordinateType> class CudaGrid;
/** \endcond */

/** \ingroup grid
    \brief Abstract wrapper of a grid.
 */
class Grid {
public:
  /** \brief Destructor */
  virtual ~Grid() {}

  /** @name Grid parameters
  @{ */

  /** \brief Dimension of the grid. */
  virtual int dim() const = 0;

  /** \brief Dimension of the space containing the grid. */
  virtual int dimWorld() const = 0;

  /** \brief Maximum level defined in this grid.

   Levels are numbered 0 ... maxLevel() with 0 the coarsest level.
   */
  virtual int maxLevel() const = 0;

  /** @}
  @name Views
  @{ */

  /** \brief View of the entities on grid level \p level. */
  virtual std::unique_ptr<GridView> levelView(size_t level) const = 0;

  /** \brief View of the leaf entities. */
  virtual std::unique_ptr<GridView> leafView() const = 0;

  /** \brief Return the topology of the grid */

  virtual GridParameters::Topology topology() const = 0;

  /** \brief Factory able to construct empty geometries of codimension 0
   compatible with this grid.

   \note For internal use.

   \todo Provide implementations for other codimensions. */
  virtual std::unique_ptr<GeometryFactory> elementGeometryFactory() const = 0;

  /** @}
  @name Others
  @{ */

  /** \brief Return a barycentrically refined grid based on the LeafView and the son map */
  virtual std::pair<shared_ptr<Grid>,Matrix<int>> barycentricGridSonPair() const = 0;

  /** \brief Return a barycentrically refined grid based on the LeafView */
  virtual shared_ptr<Grid> barycentricGrid() const = 0;

  /** \brief Return the son map for the barycentrically refined grid */
  virtual Matrix<int> barycentricSonMap() const = 0;

  /** \brief Return \p true if a barycentric refinement of this grid has
   *  been created. */
  virtual bool hasBarycentricGrid() const = 0;

//  virtual bool isBarycentricGrid() const = 0;

  /** \brief Return \p true if this grid is a barycentric representation of
   *  \p other, i.e. if this grid was created by \p other.barycentricGrid(). */
  virtual bool isBarycentricRepresentationOf(const Grid &other) const;

  /** \brief Reference to the grid's global id set. */
  virtual const IdSet &globalIdSet() const = 0;

  /** \brief Get bounding box.
   *
   *  The grid lies in the box defined by the inequalities
   *  <tt>lowerBound(i) <= x_i <= upperBound(i)</tt>, where x_i is the ith
   *  coordinate. */
  void getBoundingBox(Vector<double> &lowerBound,
                      Vector<double> &upperBound) const;

  /** \brief Get insertion index of an element. */

  virtual unsigned int elementInsertionIndex(const Entity<0> &element) const;

  /** \brief Get insertion index of a vertex for a 2d in 3d grid */

  virtual unsigned int vertexInsertionIndex(const Entity<2> &vertex) const;

  /** \brief Mark element for refinement. */

  virtual bool mark(int refCount, const Entity<0>& element) = 0;

  /** \brief Return mark status of element. */
  
  virtual int getMark(const Entity<0>& element) const = 0;

  /** \brief Pre-Adaption step */

  virtual bool preAdapt() = 0;

  /** \brief Refine mesh */

  virtual bool adapt() = 0;

  /** \brief Clean up after refinement */

  virtual void postAdapt() = 0;

  /** \brief Refine all elements refCount times. */

  virtual void globalRefine(int refCount) = 0;

  /** \brief Signal grid change to registered objects */
  void sendUpdateSignal() const;

  /** \brief Connect entity to be notified when grid updates */
  boost::signals2::connection connect(const std::function<void()>& f) const;

//  /** \brief set father of barycentric refinement */
//  virtual void setBarycentricFather(shared_ptr<Grid> fatherGrid) = 0;
//
//  /** \brief get father of barycentric refinement */
//  virtual shared_ptr<Grid> getBarycentricFather() = 0;

  /** \brief Push the mesh geometry to device memory */
  template <typename CoordinateType>
  shared_ptr<CudaGrid<CoordinateType>> pushToDevice(unsigned int deviceId) const;

private:
  /** \cond PRIVATE */
  mutable Vector<double> m_lowerBound, m_upperBound;

  mutable boost::signals2::signal<void()> gridUpdateSignal;

  mutable shared_ptr<CudaGrid<double>> cudaDoubleGridPtr;
  mutable shared_ptr<CudaGrid<float>> cudaFloatGridPtr;

  /** \endcond */
};

/** \relates Grid
 *  \brief Check whether points are inside or outside a closed grid.
 *
 *  \param[in] grid A grid representing a closed 2D surface embedded in a 3D
 *    space.
 *  \param[in] points A 2D array of dimensions (3, \c n) whose (\c i, \c j)th
 *    element is the \c i'th coordinate of \c j'th point.
 *
 *  \returns A vector of length \c n whose \c j'th element is \c true if the
 *    \c j'th point lies inside \c grid, \c false otherwise.
 *
 *  \note The results produced by this function are undefined if the grid
 *    does not represent a *closed* surface.
 *
 *  \note This function has only been tested for triangular grids.
 *
 *  \note The implementation assumes that no grid vertices are separated by less
 *    than 1e-9.
 */
std::vector<bool> areInside(const Grid &grid, const Matrix<double> &points);
std::vector<bool> areInside(const Grid &grid, const Matrix<float> &points);


} // namespace Bempp

#endif
