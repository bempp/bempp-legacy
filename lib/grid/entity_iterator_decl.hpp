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

#ifndef bempp_lib_grid_entity_iterator_decl_hpp
#define bempp_lib_grid_entity_iterator_decl_hpp

#include "entity_pointer_decl.hpp"

#include <dune/common/static_assert.hh>

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;

/** \brief Abstract base class for iterators over entities.

    \param codim Codimension of the entities iterated over

    Typical usage (ConcreteEntityIterator is some subclass of EntityIterator):

    \code
    ConcreteEntityIterator* it = ...();
    while (!it->finished())
    {
    const ConcreteEntity& e = it.entity();
    do_stuff(e);
    it->next();
    }
    delete it;
    \endcode

    \todo Various optimisations are possible in order to eliminate some virtual
    function calls. For example, next() might precompute an m_at_end flag and
    update an internal "const Entity*" pointer.

    \internal Reminder: The template parameter codim is necessary because the
    entity() method must know the codimension of the entity to which it returns
    a reference.
*/
template<int codim>
class EntityIterator: public EntityPointer<codim>
{
protected: /* Can't be changed to private, derived classes use it */
    bool m_finished;
public:
    /** \brief Increment iterator */
    virtual void next() = 0;

    /** \brief True if iterator points past the end of the iteration range */
    bool finished() const {
        return m_finished;
    }

    // virtual void reset() = 0; // Would such a method be useful?

    // This redeclaration appears so that the docstring can be changed wrt. to the base class.
    /** \brief Read-only access to the entity referenced by the iterator */
    virtual const Entity<codim>& entity() const = 0;
};

/** \brief Iterator over entities referenced by a range of Dune iterators. */
template<typename DuneEntityIt>
class ConcreteRangeEntityIterator: public EntityIterator<
    DuneEntityIt::codimension>
{
private:
    DuneEntityIt m_begin, m_end, m_cur;
    ConcreteEntity<ConcreteRangeEntityIterator::codimension,
                   typename DuneEntityIt::Entity> m_entity;

    void updateEntity() {
        if (!this->finished())
            m_entity.setDuneEntity(&*m_cur);
    }

    void updateFinished() {
        this->m_finished = (m_cur == m_end);
    }

public:
    /** Constructor. The iterator will go over the range [begin, end). */
    ConcreteRangeEntityIterator(const DuneEntityIt& begin,
                                const DuneEntityIt& end) :
        m_begin(begin), m_end(end), m_cur(begin) {
        updateFinished();
        updateEntity();
    }

    virtual void next() {
        ++m_cur;
        updateFinished();
        updateEntity();
    }

    virtual const Entity<ConcreteRangeEntityIterator::codimension>& entity() const {
        return m_entity;
    }
};

/** \brief Iterator over the subentities of codimension codim of a given entity */
template<typename DuneEntity, int codim>
class ConcreteSubentityIterator: public EntityIterator<codim>
{
    dune_static_assert(DuneEntity::codimension == 0,
                       "ConcreteSubentityIterator: only codim-0 entities "
                       "support iteration over subentities");
    dune_static_assert((int)DuneEntity::codimension < (int)ConcreteSubentityIterator::codimension,
                       "ConcreteSubentityIterator: subentity codimension "
                       "must exceed entity codimension");
public:
    typedef typename DuneEntity::template Codim<codim>::EntityPointer DuneSubentityPointer;
    typedef typename DuneSubentityPointer::Entity DuneSubentity;

private:
    const DuneEntity* m_dune_entity;
    DuneSubentityPointer m_dune_subentity_ptr;
    ConcreteEntity<ConcreteSubentityIterator::codimension, DuneSubentity> m_subentity;
    int m_cur_n;

    void updateFinished() {
        this->m_finished =
            (m_cur_n == m_dune_entity->template count<codim>());
    }

    /** \bug In order to make this compile it is necessary to define an assignment operator for
    	the nonspecialised (i.e. arbitrary codimension) Dune::FoamGridEntity class, like this:

    	(in template<int codim, int dim, class GridImp> class FoamGridEntity)
    	FoamGridEntity& operator=(const FoamGridEntity& original)
    	{
    	if (this != &original)
    	{
    	geo_.reset();
    	target_ = original.target_;
    	}
    	return *this;
    	}
    */
    void updateSubentity() {
        if (!this->finished()) {
            m_dune_subentity_ptr = m_dune_entity->template subEntity<codim>(
                m_cur_n);
            m_subentity.setDuneEntity(&*m_dune_subentity_ptr);
        }
    }

public:
    /** Constructor
    	\param dune_entity Entity whose subentities will be iterated over */
    explicit ConcreteSubentityIterator(const DuneEntity* dune_entity) :
        m_dune_entity(dune_entity), m_dune_subentity_ptr(
            m_dune_entity->template subEntity<codim>(0)), m_subentity(
                &*m_dune_subentity_ptr), m_cur_n(0) {
        updateFinished();
    }

    virtual void next() {
        ++m_cur_n;
        updateFinished();
        updateSubentity();
    }

    virtual bool finished() const {
        return m_cur_n == m_dune_entity->template count<codim>();
    }

    virtual const Entity<ConcreteSubentityIterator::codimension>& entity() const {
        return m_subentity;
    }
};

} // namespace Bempp

#endif
