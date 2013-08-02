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
#include "grid_parameters.hpp"

#include "../common/armadillo_fwd.hpp"
#include <cstddef> // size_t
#include <memory>
#include <vector>

namespace Bempp
{

/** \cond FORWARD_DECL */
template<int codim> class Entity;
class GeometryFactory;
class GridView;
class IdSet;
/** \endcond */

/** \ingroup grid
    \brief Abstract wrapper of a grid.
 */
class Grid
{
public:
    /** \brief Destructor */
    virtual ~Grid() {
    }

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
    virtual std::auto_ptr<GridView> levelView(size_t level) const = 0;

    /** \brief View of the leaf entities. */
    virtual std::auto_ptr<GridView> leafView() const = 0;

    /** \brief Return the topology of the grid */

    virtual GridParameters::Topology topology() const = 0;

    /** \brief Factory able to construct empty geometries of codimension 0
     compatible with this grid.

     \note For internal use.

     \todo Provide implementations for other codimensions. */
    virtual std::auto_ptr<GeometryFactory> elementGeometryFactory() const = 0;

    /** @}
    @name Others
    @{ */


    /** \brief Return a new barycentrically refined grid */
    virtual shared_ptr<Grid> barycentricGrid() const =0;

    /** \brief Return true of leaf-level is barycentric refinement of previous level */
    virtual bool leafIsBarycentric() const = 0;

    /** \brief Reference to the grid's global id set. */
    virtual const IdSet& globalIdSet() const = 0;

    /** \brief Get bounding box.
     *
     *  The grid lies in the box defined by the inequalities
     *  <tt>lowerBound(i) <= x_i <= upperBound(i)</tt>, where x_i is the ith
     *  coordinate. */
    void getBoundingBox(arma::Col<double>& lowerBound,
                        arma::Col<double>& upperBound) const;

private:
    /** \cond PRIVATE */
    mutable arma::Col<double> m_lowerBound, m_upperBound;
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
std::vector<bool> areInside(const Grid& grid, const arma::Mat<double>& points);
std::vector<bool> areInside(const Grid& grid, const arma::Mat<float>& points);

} // namespace Bempp

#endif
