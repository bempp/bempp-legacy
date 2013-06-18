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

#ifndef bempp_concrete_subentity_iterator_hpp
#define bempp_concrete_subentity_iterator_hpp

#include "../common/common.hpp"

#include "entity_iterator.hpp"
#include "concrete_entity_decl.hpp"

#include <dune/common/static_assert.hh>

namespace Bempp
{

/** \ingroup grid_internal
 *  \brief Iterator over the subentities of codimension \p codimSub
 *  of a given Dune entity of type \p DuneEntity. */
template<typename DuneEntity, int codimSub>
class ConcreteSubentityIterator: public EntityIterator<codimSub>
{
    dune_static_assert(DuneEntity::codimension == 0,
                       "ConcreteSubentityIterator: only codim-0 entities "
                       "support iteration over subentities");
    dune_static_assert((int)DuneEntity::codimension < (int)ConcreteSubentityIterator::codimension,
                       "ConcreteSubentityIterator: subentity codimension "
                       "must exceed entity codimension");
public:
    /** \brief Type of the Dune entity pointer corresponding to \p DuneEntity. */
    typedef typename DuneEntity::template Codim<codimSub>::EntityPointer DuneSubentityPointer;
    /** \brief Type of the appropriate subentity of \p DuneEntity. */
    typedef typename DuneSubentityPointer::Entity DuneSubentity;

private:
    const DuneEntity* m_dune_entity;
    DuneSubentityPointer m_dune_subentity_ptr;
    ConcreteEntity<ConcreteSubentityIterator::codimension, DuneSubentity> m_subentity;
    int m_cur_n;

    void updateFinished() {
        this->m_finished =
            (m_cur_n == m_dune_entity->template count<codimSub>());
    }

    void updateSubentity() {
        if (!this->finished()) {
            m_dune_subentity_ptr = m_dune_entity->template subEntity<codimSub>(
                m_cur_n);
            m_subentity.setDuneEntity(&*m_dune_subentity_ptr);
        }
    }

public:
    /** \brief Constructor.

        \param dune_entity Entity whose subentities will be iterated over */
    explicit ConcreteSubentityIterator(const DuneEntity* dune_entity) :
        m_dune_entity(dune_entity), m_dune_subentity_ptr(
            m_dune_entity->template subEntity<codimSub>(0)), m_subentity(
                &*m_dune_subentity_ptr), m_cur_n(0) {
        updateFinished();
    }

    virtual void next() {
        ++m_cur_n;
        updateFinished();
        updateSubentity();
    }

    virtual const Entity<ConcreteSubentityIterator::codimension>& entity() const {
        return m_subentity;
    }

    virtual std::auto_ptr<EntityPointer<ConcreteSubentityIterator::codimension> > frozen() const {
        const int codim = ConcreteSubentityIterator::codimension;
        return std::auto_ptr<EntityPointer<codim> >(
                    new ConcreteEntityPointer<DuneSubentityPointer>(m_subentity.duneEntity()));
    }
};

} // namespace Bempp

#endif
