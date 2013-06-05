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

#ifndef bempp_concrete_entity_decl_hpp
#define bempp_concrete_entity_decl_hpp

#include "../common/common.hpp"

#include "entity.hpp"
#include "concrete_geometry.hpp"
#include "geometry_type.hpp"

#include <boost/utility/enable_if.hpp>
#include <dune/common/static_assert.hh>
#include <dune/grid/common/geometry.hh>

namespace Bempp
{

/** \brief Wrapper of a Dune entity of type \p DuneEntity and codimension \p codim.

 \note The codimension must be given explicitly (even though it could be
 derived from the traits of \p DuneEntity) because this class needs to be
 specialised for entities of codimension 0.
 */
template<int codim, typename DuneEntity>
class ConcreteEntity: public Entity<codim>
{
    dune_static_assert((int)DuneEntity::codimension == (int)Entity<codim>::codimension,
                       "ConcreteEntity: codimension mismatch");

private:
    const DuneEntity* m_dune_entity;
    /** \internal Entity geometry. Updated on demand (on calling
     * geometry()), hence declared as mutable. */
    mutable ConcreteGeometry<typename DuneEntity::Geometry> m_geometry;

    template<typename > friend class ConcreteEntityPointer;
    template<typename, typename> friend class ConcreteRangeEntityIterator;
    template<typename, int> friend class ConcreteSubentityIterator;

    void setDuneEntity(const DuneEntity* dune_entity) {
        m_dune_entity = dune_entity;
    }

public:
    /** \brief Default constructor */
    ConcreteEntity() :
        m_dune_entity(0) {
    }

    /** \brief Constructor from a pointer to DuneEntity.

    \note This object does not acquire ownership of \p *dune_entity. */
    ConcreteEntity(const DuneEntity* dune_entity) :
        m_dune_entity(dune_entity) {
    }

    /** \brief Read-only access to the underlying Dune entity object. */
    const DuneEntity& duneEntity() const {
        return *m_dune_entity;
    }

    virtual size_t level() const {
        return m_dune_entity->level();
    }

    virtual const Geometry& geometry() const {
        if (!m_geometry.isInitialized())
            m_geometry.setDuneGeometry(m_dune_entity->geometry());
        return m_geometry;
    }

    virtual GeometryType type() const {
        return m_dune_entity->type();
    }
};

/** \brief Wrapper of a Dune entity of type \p DuneEntity and codimension 0
 */

template<typename DuneEntity>
class ConcreteEntity<0, DuneEntity> : public Entity<0>
{
    dune_static_assert((int)DuneEntity::codimension == (int)codimension,
                       "ConcreteEntity: codimension mismatch");

private:
    const DuneEntity* m_dune_entity;
    /** \internal Entity geometry. Updated on demand (on calling
     * geometry()), hence declared as mutable. */
    mutable ConcreteGeometry<typename DuneEntity::Geometry> m_geometry;

    template<typename > friend class ConcreteEntityPointer;
    template<typename, typename> friend class ConcreteRangeEntityIterator;
    template<typename, int> friend class ConcreteSubentityIterator;

    void setDuneEntity(const DuneEntity* dune_entity) {
        m_dune_entity = dune_entity;
    }

public:
    /** \brief Default constructor */
    ConcreteEntity() :
        m_dune_entity(0) {
    }

    /** \brief Constructor from a pointer to DuneEntity.

    \note This object does not acquire ownership of \p *dune_entity. */
    explicit ConcreteEntity(const DuneEntity* dune_entity) :
        m_dune_entity(dune_entity) {
    }

    /** \brief Read-only access to the underlying Dune entity object */
    const DuneEntity& duneEntity() const {
        return *m_dune_entity;
    }

    virtual size_t level() const {
        return m_dune_entity->level();
    }

    virtual const Geometry& geometry() const {
        m_geometry.setDuneGeometry(m_dune_entity->geometry());
        return m_geometry;
    }

    virtual GeometryType type() const {
        return m_dune_entity->type();
    }

    virtual std::auto_ptr<EntityPointer<0> > father() const;

    virtual bool hasFather() const {
        return m_dune_entity->hasFather();
    }

    virtual bool isLeaf() const {
        return m_dune_entity->isLeaf();
    }

    virtual bool isRegular() const {
        return m_dune_entity->isRegular();
    }

    virtual std::auto_ptr<EntityIterator<0> > sonIterator(int maxlevel) const;

    virtual bool isNew() const {
        return m_dune_entity->isNew();
    }

    virtual bool mightVanish() const {
        return m_dune_entity->mightVanish();
    }

private:
    virtual std::auto_ptr<EntityIterator<1> > subEntityCodim1Iterator() const {
        return subEntityCodimNIterator<1>();
    }
    virtual std::auto_ptr<EntityIterator<2> > subEntityCodim2Iterator() const {
        return subEntityCodimNIterator<2>();
    }
    virtual std::auto_ptr<EntityIterator<3> > subEntityCodim3Iterator() const {
        return subEntityCodimNIterator<3>();
    }

    // these methods are implemented in entity.hpp (outside the declaration of
    // ConcreteEntity) because they need to know the full declaration of
    // concrete iterator (which may not be available at this stage)
    template <int codimSub>
    typename boost::disable_if_c<codimSub <= DuneEntity::dimension, std::auto_ptr<EntityIterator<codimSub> > >::type
    subEntityCodimNIterator() const;

    template <int codimSub>
    typename boost::enable_if_c<codimSub <= DuneEntity::dimension, std::auto_ptr<EntityIterator<codimSub> > >::type
    subEntityCodimNIterator() const;

    virtual size_t subEntityCodim1Count() const {
        return subEntityCodimNCount<1>();
    }
    virtual size_t subEntityCodim2Count() const {
        return subEntityCodimNCount<2>();
    }
    virtual size_t subEntityCodim3Count() const {
        return subEntityCodimNCount<3>();
    }

    template <int codimSub>
    typename boost::disable_if_c<codimSub <= DuneEntity::dimension, size_t>::type
    subEntityCodimNCount() const {
        return 0;
    }

    template <int codimSub>
    typename boost::enable_if_c<codimSub <= DuneEntity::dimension, size_t>::type
    subEntityCodimNCount() const {
        return m_dune_entity->template count<codimSub>();
    }
};

} // namespace Bempp

#endif
