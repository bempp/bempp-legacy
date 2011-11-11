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

#ifndef bempp_lib_grid_3d_entity_hpp
#define bempp_lib_grid_3d_entity_hpp

#include "entity_decl.hpp"
#include "entity_iterator.hpp"
#include "entity_pointer.hpp"

namespace Bempp {
namespace ThreeD {

template<typename DuneEntity>
inline EntityIterator<1>* ConcreteEntity<0, DuneEntity>::subEdgeIterator() const {
	const int codim = 1;
	typedef ConcreteSubentityIterator<DuneEntity, codim> ConcIterator;
	return new ConcIterator(m_dune_entity);
}

template<typename DuneEntity>
inline EntityIterator<2>* ConcreteEntity<0, DuneEntity>::subVertexIterator() const {
	const int codim = 2;
	typedef ConcreteSubentityIterator<DuneEntity, codim> ConcIterator;
	return new ConcIterator(m_dune_entity);
}

//  virtual SurfaceGridIntersectionIterator ileafbegin() const
//  {
//    typedef typename GridType::LeafIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ileafbegin();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

//  virtual SurfaceGridIntersectionIterator ileafend() const
//  {
//    typedef typename GridType::LeafIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ileafend();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

//  virtual SurfaceGridIntersectionIterator ilevelbegin() const
//  {
//    typedef typename GridType::LevelIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ilevelbegin();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

//  virtual SurfaceGridIntersectionIterator ilevelend() const
//  {
//    typedef typename GridType::LevelIntersectionIterator DuneIterator;
//    typedef ConcreteSurfaceGridIntersectionIterator<DuneIterator> ConcreteIterator;
//    DuneIterator dit = m_dune_entity->ilevelend();
//    return SurfaceGridIntersectionIterator(new ConcreteIterator(dit));
//  }

template<typename DuneEntity>
inline EntityPointer<0>* ConcreteEntity<0, DuneEntity>::father() const {
	typedef ConcreteEntityPointer<typename DuneEntity::EntityPointer> ConcPointer;
	return new ConcPointer(m_dune_entity->father());
}

template<typename DuneEntity>
inline EntityIterator<0>* ConcreteEntity<0, DuneEntity>::sonIterator(
		int maxlevel) const {
	typedef typename DuneEntity::HierarchicIterator DuneIterator;
	typedef ConcreteRangeEntityIterator<DuneIterator> ConcIterator;
	return new ConcIterator(m_dune_entity->hbegin(maxlevel),
			m_dune_entity->hend(maxlevel));
}

} // namespace Bempp
} // namespace ThreeD

#endif
