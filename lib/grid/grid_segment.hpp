// Copyright (C) 2011-2013 by the BEM++ Authors
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

#ifndef bempp_grid_segment_hpp
#define bempp_grid_segment_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include <boost/array.hpp>
#include <set>
#include <vector>

namespace Bempp {

class Grid;

/** \brief Segment of a grid.
 *
 *  This class represent a segment (part) of a grid. You can use it to create
 *  spaces of functions defined only on subsets of grids, for example to impose
 *  different boundary conditions on different parts of a single grid.
 *
 *  Technically, a grid segment is defined by means of the sets of indices of
 *  entities of different codimensions (vertices, edges and elements, for 2D
 *  grids) that it contains.
 */
class GridSegment {
public:
  /** \brief Return a GridSegment representing the whole grid \p grid at
      a given level (level = -1 for leafView). */
  static GridSegment wholeGrid(const Grid &grid, int level = 0);

  /** \brief Return a GridSegment representing an open domain of a grid.
   *
   *  This function returns a GridSegment representing the domain with index
   *  \p domain of the grid \p grid, not including its boundary. The domain
   *  index corresponds to the <it>physical entity index</it> in Gmsh files.*/
  static GridSegment openDomain(const Grid &grid, int domain, int level = 0);

  /** \brief Overload to specify multiple domains. */
  static GridSegment openDomain(const Grid &grid, 
          const std::vector<int>& domains, int level = 0);

  /** \brief Return a GridSegment representing a closed domain of a grid.
   *
   *  This function returns a GridSegment representing the domain with index
   *  \p domain of the grid \p grid, including its boundary. The domain
   *  index corresponds to the <it>physical entity index</it> in Gmsh files.*/
  static GridSegment closedDomain(const Grid &grid, int domain, int level = 0);

  /** \brief Overload to specify multiple domains. */
  static GridSegment closedDomain(const Grid &grid, 
          const std::vector<int>& domains, int level = 0);

  /** \brief Constructor.
   *
   *  Construct a new segment of the grid \p grid including all its
   *  constituent entities except those listed in \p excludedEntitiesCodim0,
   *  \p excludedEntitiesCodim1 and so on (each set contains the indices of
   *  entities of a particular codimension that should be excluded from the
   *  segment). The indices should correspond to those defined by the leaf
   *  view of the grid. The \p grid reference is used only in the constructor
   *  and is not stored in the constructed GridSegment.
   */
  GridSegment(const Grid &grid, const std::set<int> &excludedEntitiesCodim0,
              const std::set<int> &excludedEntitiesCodim1,
              const std::set<int> &excludedEntitiesCodim2,
              const std::set<int> &excludedEntitiesCodim3, int level = 0);

  /** \brief Constructor.
   *
   *  Construct a new segment of a grid including all its constituent
   *  entities except those listed in \p excludedEntitiesCodim0, \p
   *  excludedEntitiesCodim1 and so on (each set contains the indices of
   *  entities of a particular codimension that should be excluded from the
   *  segment). The indices should correspond to those defined by the leaf
   *  view of the grid. The arguments \p entityCountCodim0, \p
   *  entityCountCodim1 and so on should contain the number of entities of a
   *  particular codimension contained in the grid. */
  GridSegment(int entityCountCodim0, int entityCountCodim1,
              int entityCountCodim2, int entityCountCodim3,
              const std::set<int> &excludedEntitiesCodim0,
              const std::set<int> &excludedEntitiesCodim1,
              const std::set<int> &excludedEntitiesCodim2,
              const std::set<int> &excludedEntitiesCodim3);

  /** \brief Return true if the segment contains the entity of codimension \p
  codim and index \p index, false otherwise. */
  bool contains(int codim, int index) const;

  /** \brief Mark array elements corresponding to entities that do not belong
   *  to this GridSegment.
   *
   *  This function assigns the value \p mark to all the elements of the
   *  vector \p marks whose indices correspond to entities of codimension \p
   *  codim that do not belong to this GridSegment. All other elements are
   *  set to 0.
   *
   *  \note This function may be removed in future versions of BEM++. */
  void markExcludedEntities(int codim, std::vector<int> &marks,
                            int mark = -1) const;

  /** \brief Return the complement of this grid segment.
   *
   *  The returned GridSegment contains those and only those entities that
   *  do not belong to this GridSegment. */
  GridSegment complement() const;

  /** \brief Return the union of this GridSegment and the GridSegment
   *  \p other.
   *
   *  The returned GridSegment contains those and only those entities that
   *  belong to this GridSegment and/or the GridSegment \p other. */
  GridSegment union_(const GridSegment &other) const;

  /** \brief Return the difference of this GridSegment and the GridSegment
   *  \p other.
   *
   *  The returned GridSegment contains those and only those entities that
   *  belong to this GridSegment, but do not belong to the GridSegment \p
   *  other. */
  GridSegment difference(const GridSegment &other) const;

  /** \brief Return the intersection of this GridSegment and the GridSegment
   *  \p other.
   *
   *  The returned GridSegment contains those and only those entities that
   *  belong to both this GridSegment and the GridSegment \p other. */
  GridSegment intersection(const GridSegment &other) const;

private:
  boost::array<int, 4> m_entityCounts;
  boost::array<std::set<int>, 4> m_excludedEntities;
};

GridSegment gridSegmentWithPositiveX(const Grid &grid, int level = 0);

} // namespace Bempp

#endif
