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

#ifndef bempp_concrete_index_set_hpp
#define bempp_concrete_index_set_hpp

#include "index_set.hpp"
#include "concrete_entity.hpp"

namespace Bempp
{


/** \brief Wrapper of the index set specific to a Dune grid view class \p DuneGridView

 \internal The grid view class, rather than an index set class, is used as a
 template parameter because the latter doesn't provide information about the
 entity type.

 For consistency with \p IdSet it would be possible to take as parameters
 \p DuneGrid and \p DuneIndexSet instead.
 */
template<typename DuneGridView>
class ConcreteIndexSet: public IndexSet
{
public:
    /** \brief Type of the wrapped Dune index set. */
    typedef typename DuneGridView::IndexSet DuneIndexSet;

private:
    const DuneIndexSet* m_dune_index_set;

public:
    /** \brief Constructor.

    This object does not assume ownership of \p *dune_index_set.*/
    explicit ConcreteIndexSet(const DuneIndexSet* dune_index_set) :
        m_dune_index_set(dune_index_set) {
    }

    /** \brief Read-only access to the underlying Dune index set. */
    const DuneIndexSet& duneIndexSet() const {
        return *m_dune_index_set;
    }

    virtual IndexType entityIndex(const Entity<0>& e) const {
        return entityCodimNIndex(e);
    }
    virtual IndexType entityIndex(const Entity<1>& e) const {
        return entityCodimNIndex(e);
    }
    virtual IndexType entityIndex(const Entity<2>& e) const {
        return entityCodimNIndex(e);
    }
    virtual IndexType entityIndex(const Entity<3>& e) const {
        return entityCodimNIndex(e);
    }

private:
    template <int codim>
    typename boost::disable_if_c<codim <= DuneGridView::dimension, IndexType>::type
    entityCodimNIndex(const Entity<codim>& e) const {
        throw std::logic_error("IndexSet::entityIndex(): invalid entity codimension");
    }

    template <int codim>
    typename boost::enable_if_c<codim <= DuneGridView::dimension, IndexType>::type
    entityCodimNIndex(const Entity<codim>& e) const {
        typedef typename DuneGridView::template Codim<codim>::Entity DuneEntity;
        typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_index_set->index(ce.duneEntity());
    }
};

} // namespace Bempp

#endif
