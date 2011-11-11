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

#ifndef bempp_lib_grid_3d_entity_pointer_decl_hpp
#define bempp_lib_grid_3d_entity_pointer_decl_hpp

namespace Bempp {

// Forward declarations
template<int codim> class Entity;

/**
 \brief Abstract base class for an object providing read-only access to an
 entity of codimension codim.

 \internal Reminder: The template parameter codim is necessary because the
 entity() method must know the codimension of the entity to which it returns
 a reference.
 */
template<int codim>
class EntityPointer {
public:
	/** Destructor */
	virtual ~EntityPointer() {
	}

	/** Entity codimension */
	enum {
		codimension = codim
	};

	/** Read-only access to the underlying entity */
	virtual const Entity<codim>& entity() const = 0;
};

/**
 \brief A class providing read-only access to a wrapper of a Dune entity of
 type DuneEntity.
 */
template<typename DuneEntityPointer>
class ConcreteEntityPointer: public EntityPointer<DuneEntityPointer::codimension> {
protected:
	typedef typename DuneEntityPointer::Entity DuneEntity;
	DuneEntityPointer m_dune_entity_ptr;
	ConcreteEntity<ConcreteEntityPointer::codimension, DuneEntity> m_entity;

	void updateEntity() {
		m_entity.setDuneEntity(&*m_dune_entity_ptr);
	}

public:
	/** Constructor */
	explicit ConcreteEntityPointer(const DuneEntityPointer& dune_entity_pointer) :
			m_dune_entity_ptr(dune_entity_pointer) {
		updateEntity();
	}

	virtual const Entity<DuneEntityPointer::codimension>& entity() const {
		return m_entity;
	}
};

} // namespace Bempp

#endif
