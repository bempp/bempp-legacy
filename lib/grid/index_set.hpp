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

#ifndef bempp_lib_grid_3d_index_set_hpp
#define bempp_lib_grid_3d_index_set_hpp

#include "index_set_decl.hpp"
#include "entity.hpp"

namespace Bempp {

template<typename DuneGridView>
inline IndexSet::IndexType ConcreteIndexSet<DuneGridView>::faceIndex(
		const Entity<0>& entity) const {
	const int codim = 0;
	typedef typename DuneGridView::template Codim<codim>::Entity DuneEntity;
	typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
	const ConcEntity& ce = dynamic_cast<const ConcEntity&>(entity);
	return m_dune_index_set->index(ce.duneEntity());
}

template<typename DuneGridView>
inline IndexSet::IndexType ConcreteIndexSet<DuneGridView>::edgeIndex(
		const Entity<1>& entity) const {
	const int codim = 1;
	typedef typename DuneGridView::template Codim<codim>::Entity DuneEntity;
	typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
	const ConcEntity& ce = dynamic_cast<const ConcEntity&>(entity);
	return m_dune_index_set->index(ce.duneEntity());
}

template<typename DuneGridView>
inline IndexSet::IndexType ConcreteIndexSet<DuneGridView>::vertexIndex(
		const Entity<2>& entity) const {
	const int codim = 2;
	typedef typename DuneGridView::template Codim<codim>::Entity DuneEntity;
	typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
	const ConcEntity& ce = dynamic_cast<const ConcEntity&>(entity);
	return m_dune_index_set->index(ce.duneEntity());
}

} // namespace Bempp

#endif
