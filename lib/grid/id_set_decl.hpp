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

#ifndef bempp_id_set_decl_hpp
#define bempp_id_set_decl_hpp

#include <boost/utility/enable_if.hpp>

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;

/** \brief Abstract wrapper of an id set. */
class IdSet
{
public:
    /** \brief Destructor. */
    virtual ~IdSet() {
    }

    /** \brief Id type.

     \internal Sadly, it is necessary to specify this type uniformly for all grid classes.
     */
    typedef unsigned int IdType;

    /** \brief Id of the entity \p e of codimension 0. */
    virtual IdType entityId(const Entity<0>& e) const = 0;
    /** \brief Id of the entity \p e of codimension 1. */
    virtual IdType entityId(const Entity<1>& e) const = 0;
    /** \brief Id of the entity \p e of codimension 2. */
    virtual IdType entityId(const Entity<2>& e) const = 0;
    /** \brief Id of the entity \p e of codimension 3. */
    virtual IdType entityId(const Entity<3>& e) const = 0;
};

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
    entityCodimNId(const Entity<codim>& e) const;
    template <int codim>
    typename boost::enable_if_c<codim <= DuneGrid::dimension, IdType>::type
    entityCodimNId(const Entity<codim>& e) const;
};

} // namespace Bempp

#endif // ID_SET_HPP
