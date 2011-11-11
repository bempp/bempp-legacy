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

#ifndef bempp_lib_grid_3d_id_set_decl_hpp
#define bempp_lib_grid_3d_id_set_decl_hpp

namespace Bempp {
namespace ThreeD {

// Forward declarations
template<int codim> class Entity;

/** Abstract wrapper of an id set */
class IdSet {
public:
	/** Detructor */
	virtual ~IdSet() {
	}

	/** Id type

	 \internal Sadly, it is necessary to specify this type uniformly for all grid classes. 
	 */
	typedef unsigned int IdType;

	/** \brief Map face entity to id. */
	virtual IdType faceId(const Entity<0>& entity) const = 0;

	/** \brief Map edge entity to id. */
	virtual IdType edgeId(const Entity<1>& entity) const = 0;

	/** \brief Map vertex entity to id. */
	virtual IdType vertexId(const Entity<2>& entity) const = 0;
};

/** \brief Wrapper of the Dune id set of type DuneIdSet providing access to the
 entities of a Dune grid of type DuneGrid.

 \internal Both these typenames are needed because Dune::IdSet does not export
 entity type.
 */
template<typename DuneGrid, typename DuneIdSet>
class ConcreteIdSet: public IdSet {
protected:
	const DuneIdSet* m_dune_id_set;

public:
	/**  Constructor */
	explicit ConcreteIdSet(const DuneIdSet* dune_id_set) :
			m_dune_id_set(dune_id_set) {
	}

	/** Read-only access to the underlying Dune id set */
	const DuneIdSet& duneIdSet() const {
		return *m_dune_id_set;
	}

	virtual IdType faceId(const Entity<0>& entity) const;
	virtual IdType edgeId(const Entity<1>& entity) const;
	virtual IdType vertexId(const Entity<2>& entity) const;
};

} // namespace Bempp
} // namespace ThreeD

#endif // ID_SET_HPP
