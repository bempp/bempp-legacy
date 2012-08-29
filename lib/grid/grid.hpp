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

#include "../common/common.hpp"
#include "grid_parameters.hpp"

#include <memory>
#include <cstddef> // size_t

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;
class GeometryFactory;
class GridView;
class IdSet;

/** \brief Abstract wrapper of a grid.

    Functions related to parallelization are not wrapped yet.
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

    /** \brief Number of boundary segments within the macro (level-0) grid. */
    virtual size_t boundarySegmentCount() const = 0;

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
    @name Id sets
    @{ */

    /** \brief Reference to the grid's global id set. */
    virtual const IdSet& globalIdSet() const = 0;
};

} // namespace Bempp

#endif
