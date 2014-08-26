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

#include "grid_segment.hpp"

#include "../common/acc.hpp"

#include "entity.hpp"
#include "entity_iterator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "grid_view.hpp"
#include "index_set.hpp"

#include <stdexcept>

namespace Bempp {

namespace {

template <int codim>
std::set<int> entitiesWithNonpositiveX(const GridView &view) {
  std::set<int> result;
  std::unique_ptr<EntityIterator<codim>> it = view.entityIterator<codim>();
  const IndexSet &indexSet = view.indexSet();
  arma::Col<double> center;
  while (!it->finished()) {
    const Entity<codim> &entity = it->entity();
    entity.geometry().getCenter(center);
    if (center(0) <= 0.)
      result.insert(indexSet.entityIndex(entity));
    it->next();
  }
  return result;
}

} // namespace

GridSegment::GridSegment(const Grid &grid,
                         const std::set<int> &excludedEntitiesCodim0,
                         const std::set<int> &excludedEntitiesCodim1,
                         const std::set<int> &excludedEntitiesCodim2,
                         const std::set<int> &excludedEntitiesCodim3,
                         int level) {

  std::unique_ptr<GridView> view;
  if (level == -1)
    view = grid.leafView();
  else
    view = grid.levelView(level);
  for (int i = 0; i < 4; ++i)
    m_entityCounts[i] = view->entityCount(i);
  m_excludedEntities[0] = excludedEntitiesCodim0;
  m_excludedEntities[1] = excludedEntitiesCodim1;
  m_excludedEntities[2] = excludedEntitiesCodim2;
  m_excludedEntities[3] = excludedEntitiesCodim3;
}

GridSegment::GridSegment(int entityCountCodim0, int entityCountCodim1,
                         int entityCountCodim2, int entityCountCodim3,
                         const std::set<int> &excludedEntitiesCodim0,
                         const std::set<int> &excludedEntitiesCodim1,
                         const std::set<int> &excludedEntitiesCodim2,
                         const std::set<int> &excludedEntitiesCodim3) {
  m_entityCounts[0] = entityCountCodim0;
  m_entityCounts[1] = entityCountCodim1;
  m_entityCounts[2] = entityCountCodim2;
  m_entityCounts[3] = entityCountCodim3;
  m_excludedEntities[0] = excludedEntitiesCodim0;
  m_excludedEntities[1] = excludedEntitiesCodim1;
  m_excludedEntities[2] = excludedEntitiesCodim2;
  m_excludedEntities[3] = excludedEntitiesCodim3;
}

GridSegment GridSegment::wholeGrid(const Grid &grid, int level) {
  std::set<int> emptySet;
  return GridSegment(grid, emptySet, emptySet, emptySet, emptySet, level);
}

GridSegment GridSegment::openDomain(const Grid &grid, int domain, int level) {
  const int gridDim = grid.dim();
  std::unique_ptr<GridView> view;
  if (level == -1)
    view = grid.leafView();
  else
    view = grid.levelView(level);
  const IndexSet &indexSet = view->indexSet();

  boost::array<std::vector<bool>, 4> entirelyInDomain;
  for (int codim = 1; codim <= gridDim; ++codim)
    entirelyInDomain[codim].resize(view->entityCount(codim), true);

  boost::array<std::set<int>, 4> excludedEntities;
  std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    if (e.domain() != domain) {
      excludedEntities[0].insert(indexSet.entityIndex(e));
      for (int codim = 1; codim <= gridDim; ++codim) {
        int subEntityCount =
            (codim == 1) ? e.subEntityCount<1>() : (codim == 2)
                                                       ? e.subEntityCount<2>()
                                                       : e.subEntityCount<3>();
        for (int i = 0; i < subEntityCount; ++i) {
          int index = indexSet.subEntityIndex(e, i, codim);
          acc(entirelyInDomain[codim], index) = false;
        }
      }
    }
    it->next();
  }
  for (int codim = 1; codim <= gridDim; ++codim)
    for (int index = 0; index < entirelyInDomain[codim].size(); ++index)
      if (!acc(entirelyInDomain[codim], index))
        excludedEntities[codim].insert(index);
  return GridSegment(view->entityCount(0), view->entityCount(1),
                     view->entityCount(2), view->entityCount(3),
                     excludedEntities[0], excludedEntities[1],
                     excludedEntities[2], excludedEntities[3]);
}

GridSegment GridSegment::closedDomain(const Grid &grid, int domain, int level) {
  const int gridDim = grid.dim();
  std::unique_ptr<GridView> view;
  if (level == -1)
    view = grid.leafView();
  else
    view = grid.levelView(level);
  const IndexSet &indexSet = view->indexSet();

  boost::array<std::vector<bool>, 4> adjacentToDomain;
  for (int codim = 1; codim <= gridDim; ++codim)
    adjacentToDomain[codim].resize(view->entityCount(codim), false);

  boost::array<std::set<int>, 4> excludedEntities;
  std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    if (e.domain() != domain)
      excludedEntities[0].insert(indexSet.entityIndex(e));
    else {
      for (int codim = 1; codim <= gridDim; ++codim) {
        int subEntityCount =
            (codim == 1) ? e.subEntityCount<1>() : (codim == 2)
                                                       ? e.subEntityCount<2>()
                                                       : e.subEntityCount<3>();
        for (int i = 0; i < subEntityCount; ++i) {
          int index = indexSet.subEntityIndex(e, i, codim);
          acc(adjacentToDomain[codim], index) = true;
        }
      }
    }
    it->next();
  }
  for (int codim = 1; codim <= gridDim; ++codim)
    for (int index = 0; index < adjacentToDomain[codim].size(); ++index)
      if (!acc(adjacentToDomain[codim], index))
        excludedEntities[codim].insert(index);
  return GridSegment(view->entityCount(0), view->entityCount(1),
                     view->entityCount(2), view->entityCount(3),
                     excludedEntities[0], excludedEntities[1],
                     excludedEntities[2], excludedEntities[3]);
}

bool GridSegment::contains(int codim, int index) const {
  if (codim < 0 || codim > 3)
    throw std::invalid_argument("GridSegment::contains(): codim must be "
                                "0, 1, 2 or 3");
  return index >= 0 && index < m_entityCounts[codim] &&
         m_excludedEntities[codim].find(index) ==
             m_excludedEntities[codim].end();
}

void GridSegment::markExcludedEntities(int codim, std::vector<int> &marks,
                                       int mark) const {
  if (codim < 0 || codim > 3)
    throw std::invalid_argument("GridSegment::begin(): codim must be "
                                "0, 1, 2 or 3");
  marks.resize(m_entityCounts[codim]);
  std::fill(marks.begin(), marks.end(), 0);
  for (std::set<int>::const_iterator it = m_excludedEntities[codim].begin();
       it != m_excludedEntities[codim].end(); ++it)
    acc(marks, *it) = mark;
}

GridSegment GridSegment::complement() const {
  boost::array<std::set<int>, 4> includedEntities;
  for (int codim = 0; codim < 4; ++codim) {
    std::set<int>::const_iterator it = m_excludedEntities[codim].begin();
    for (int index = 0; index < m_entityCounts[codim]; ++index) {
      if (it == m_excludedEntities[codim].end() || index < *it)
        includedEntities[codim].insert(index);
      else
        ++it;
    }
  }
  return GridSegment(m_entityCounts[0], m_entityCounts[1], m_entityCounts[2],
                     m_entityCounts[3], includedEntities[0],
                     includedEntities[1], includedEntities[2],
                     includedEntities[3]);
}

GridSegment GridSegment::union_(const GridSegment &other) const {
  boost::array<std::set<int>, 4> excludedEntities;
  for (int codim = 0; codim < 4; ++codim)
    std::set_intersection(m_excludedEntities[codim].begin(),
                          m_excludedEntities[codim].end(),
                          other.m_excludedEntities[codim].begin(),
                          other.m_excludedEntities[codim].end(),
                          std::inserter(excludedEntities[codim],
                                        excludedEntities[codim].begin()));
  return GridSegment(m_entityCounts[0], m_entityCounts[1], m_entityCounts[2],
                     m_entityCounts[3], excludedEntities[0],
                     excludedEntities[1], excludedEntities[2],
                     excludedEntities[3]);
}

GridSegment GridSegment::difference(const GridSegment &other) const {
  GridSegment otherC = other.complement();
  return intersection(otherC);
}

GridSegment GridSegment::intersection(const GridSegment &other) const {
  boost::array<std::set<int>, 4> excludedEntities;
  for (int codim = 0; codim < 4; ++codim)
    std::set_union(m_excludedEntities[codim].begin(),
                   m_excludedEntities[codim].end(),
                   other.m_excludedEntities[codim].begin(),
                   other.m_excludedEntities[codim].end(),
                   std::inserter(excludedEntities[codim],
                                 excludedEntities[codim].begin()));
  return GridSegment(m_entityCounts[0], m_entityCounts[1], m_entityCounts[2],
                     m_entityCounts[3], excludedEntities[0],
                     excludedEntities[1], excludedEntities[2],
                     excludedEntities[3]);
}

GridSegment gridSegmentWithPositiveX(const Grid &grid, int level) {
  boost::array<std::set<int>, 4> excludedEntities;
  std::unique_ptr<GridView> view;
  if (level == -1)
    view = grid.leafView();
  else
    view = grid.levelView(level);
  excludedEntities[0] = entitiesWithNonpositiveX<0>(*view);
  excludedEntities[1] = entitiesWithNonpositiveX<1>(*view);
  excludedEntities[2] = entitiesWithNonpositiveX<2>(*view);

  return GridSegment(grid, excludedEntities[0], excludedEntities[1],
                     excludedEntities[2], excludedEntities[3]);
}

} // namespace Bempp
