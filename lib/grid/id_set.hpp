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

#ifndef bempp_lib_grid_3d_id_set_hpp
#define bempp_lib_grid_3d_id_set_hpp

#include "id_set_decl.hpp"
#include "entity.hpp"

namespace Bempp
{
namespace ThreeD
{

template<typename DuneGrid, typename DuneIdSet>
inline IdSet::IdType ConcreteIdSet<DuneGrid, DuneIdSet>::faceId(const Entity<0>& entity) const
{
  // errors here -> ConcreteIdSet needs to be passed both DuneGrid and IdSet (bc. DuneGrid has two id sets -- Global and Local -- but those don't provide info abouyt entity type)
  // for a start maybe it's better to comment out these things.
  const int codim = 0;
  typedef typename DuneGrid::template Codim<codim>::Entity DuneEntity;
  typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
  const ConcEntity& ce = dynamic_cast<const ConcEntity&>(entity);
  return m_dune_id_set->id(ce.duneEntity());
  return 0;
}

template<typename DuneGrid, typename DuneIdSet>
inline IdSet::IdType ConcreteIdSet<DuneGrid, DuneIdSet>::edgeId(const Entity<1>& entity) const
{
  const int codim = 1;
  typedef typename DuneGrid::template Codim<codim>::Entity DuneEntity;
  typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
  const ConcEntity& ce = dynamic_cast<const ConcEntity&>(entity);
  return m_dune_id_set->id(ce.duneEntity());
  return 0;
}

template<typename DuneGrid, typename DuneIdSet>
inline IdSet::IdType ConcreteIdSet<DuneGrid, DuneIdSet>::vertexId(const Entity<2>& entity) const
{
  const int codim = 2;
  typedef typename DuneGrid::template Codim<codim>::Entity DuneEntity;
  typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
  const ConcEntity& ce = dynamic_cast<const ConcEntity&>(entity);
  return m_dune_id_set->id(ce.duneEntity());
  return 0;
}

} // namespace Bempp
} // namespace ThreeD

#endif // ID_SET_HPP
