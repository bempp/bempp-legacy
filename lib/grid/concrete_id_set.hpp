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

#ifndef bempp_concrete_id_set_hpp
#define bempp_concrete_id_set_hpp

#include "id_set.hpp"
#include "concrete_entity.hpp"

#include <stdexcept>

namespace Bempp
{

/** \brief Wrapper of a Dune id set of type \p DuneIdSet providing access to the
 entities of a Dune grid of type \p DuneGrid.

 \internal Both these typenames are needed because <tt>Dune::IdSet</tt> does not export
 entity type.
 */
template<typename DuneGrid, typename DuneIdSet>
class ConcreteIdSet: public IdSet
{
private:
    const DuneIdSet* m_dune_id_set;

public:
    /** \brief Constructor.

      This object does not assume ownership of \p *dune_id_set.
    */
    explicit ConcreteIdSet(const DuneIdSet* dune_id_set) :
        m_dune_id_set(dune_id_set) {
    }

    /** \brief Read-only access to the underlying Dune id set. */
    const DuneIdSet& duneIdSet() const {
        return *m_dune_id_set;
    }

    virtual IdType entityId(const Entity<0>& e) const {
        return entityCodimNId(e);
    }
    virtual IdType entityId(const Entity<1>& e) const {
        return entityCodimNId(e);
    }
    virtual IdType entityId(const Entity<2>& e) const {
        return entityCodimNId(e);
    }
    virtual IdType entityId(const Entity<3>& e) const {
        return entityCodimNId(e);
    }

    virtual IdType subEntityId(const Entity<0>& e, int i, unsigned int codimSub) const {
#ifndef NDEBUG
        // Prevent an assert in FoamGrid from crashing the Python interpreter
        if (codimSub > DuneGrid::dimension)
            throw std::invalid_argument("IndexSet::subEntityIndex(): codimSub exceeds grid dimension");
#endif
        typedef typename DuneGrid::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_id_set->subId(ce.duneEntity(), i, codimSub);
    }

private:
    // Below there are essentially two overloads of
    // template <int codim>
    // IdType entityCodimNId(const Entity<codim>& e) const;
    // valid for codim > DuneGrid::dimension (throws exception)
    // and for the opposite case (does real work).
    //
    // This is not a very pretty code, ideas how to improve it are welcome!
    // The problem is that we cannot allow the compiler to instantiate
    // Dune::Entity objects with codim > DuneGrid::dimension
    // because that throws a static assert in Dune code.
    template <int codim>
    typename boost::disable_if_c<codim <= DuneGrid::dimension, IdType>::type
    entityCodimNId(const Entity<codim>& e) const {
        throw std::logic_error("IdSet::entityId(): invalid entity codimension");
    }

    template <int codim>
    typename boost::enable_if_c<codim <= DuneGrid::dimension, IdType>::type
    entityCodimNId(const Entity<codim>& e) const {
        typedef typename DuneGrid::template Codim<codim>::Entity DuneEntity;
        typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_id_set->id(ce.duneEntity());
    }
};

} // namespace Bempp

#endif // ID_SET_HPP
