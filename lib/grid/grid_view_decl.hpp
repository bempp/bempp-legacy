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

#ifndef bempp_lib_grid_3d_grid_view_decl_hpp
#define bempp_lib_grid_3d_grid_view_decl_hpp

#include "index_set_decl.hpp"

namespace Dune {
class GeometryType;
}

namespace Bempp {

// Forward declarations
template<int codim> class Entity;
template<int codim> class EntityIterator;
class IndexSet;

/** Abstract wrapper of a grid view */
class GridView {
public:
	/** Destructor */
	virtual ~GridView() {
	}

	// TODO: return base grid
	// virtual & grid() const = 0;

	/** \brief The index set */
	virtual const IndexSet& indexSet() const = 0;

	/** \brief Number of entities in a given codimension */
	virtual int size(int codim) const = 0;

	/** \brief Number of entities with a given geometry type */
	virtual int size(const Dune::GeometryType &type) const = 0;

	/** \brief True if the given vertex is contained in this grid view.
	 *
	 * \note If e is not an element of the grid, then
	 *       the result of containsVertex() is undefined.
	 */
	virtual bool containsVertex(const Entity<2>& e) const = 0;

	/** \brief True if the given edge is contained in this grid view.
	 *
	 * \note If e is not an element of the grid, then
	 *       the result of containsEdge() is undefined.
	 */
	virtual bool containsEdge(const Entity<1>& e) const = 0;

	/** \brief True if the given face is contained in this grid view.
	 *
	 * \note If e is not an element of the grid, then
	 *       the result of containsFace() is undefined.
	 */
	virtual bool containsFace(const Entity<0>& e) const = 0;

	/** \brief Iterator over faces contained in this view.
	 *
	 *  The caller is responsible for freeing the returned pointer.
	 */
	virtual EntityIterator<0>* faceIterator() const = 0;

	/** \brief Iterator over edges contained in this view.
	 *
	 *  The caller is responsible for freeing the returned pointer.
	 */
	virtual EntityIterator<1>* edgeIterator() const = 0;

	/** \brief Iterator over vertices contained in this view.
	 *
	 *  The caller is responsible for freeing the returned pointer.
	 */
	virtual EntityIterator<2>* vertexIterator() const = 0;

	// TODO: Intersection iterators.
//  virtual Iterator
//  ibegin(const Entity<0>& entity) const = 0;
//  virtual Iterator
//  iend(const Entity<0>& entity) const = 0;
};

/** Wrapper of a Dune grid view of type DuneGridView. */
template<typename DuneGridView>
class ConcreteGridView: public GridView {
protected:
	DuneGridView m_dune_gv;
	ConcreteIndexSet<DuneGridView> m_index_set;

public:
	/** Constructor */
	explicit ConcreteGridView(const DuneGridView& dune_gv) :
			m_dune_gv(dune_gv), m_index_set(&dune_gv.indexSet()) {
	}

	/** Access to the underlying Dune grid view object. Use at your own risk! */
	const DuneGridView& duneGridView() const {
		return m_dune_gv;
	}
	/** Read-only access to the underlying Dune grid view object. */
	DuneGridView& duneGridView() {
		return m_dune_gv;
	}

	// virtual & grid() const
	//   { return Concrete<const DuneGridView::Grid>(&m_dune_gv.grid()); }

	virtual const IndexSet& indexSet() const {
		return m_index_set;
	}

	virtual int size(int codim) const {
		return m_dune_gv.size(codim);
	}

	virtual int size(const Dune::GeometryType &type) const {
		return m_dune_gv.size(type);
	}

	virtual bool containsFace(const Entity<0>& e) const;
	virtual bool containsEdge(const Entity<1>& e) const;
	virtual bool containsVertex(const Entity<2>& e) const;
	virtual EntityIterator<0>* faceIterator() const;
	virtual EntityIterator<1>* edgeIterator() const;
	virtual EntityIterator<2>* vertexIterator() const;
};

} // namespace Bempp

#endif
