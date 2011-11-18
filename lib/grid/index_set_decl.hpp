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

#ifndef bempp_index_set_decl_hpp
#define bempp_index_set_decl_hpp

namespace Bempp
{

// Forward declarations
template<int codim> class Entity;

/** Abstract wrapper of an index set */
class IndexSet
{
public:
    /** Detructor */
    virtual ~IndexSet() {
    }

    /** Index type

     \internal Sadly, it is necessary to specify this type uniformly for all grid classes.
     */
    typedef unsigned int IndexType;

    /** \brief Map entity of codimension 0 to index.

     The result of calling this method with an entity that is not
     in the index set is undefined.

     \return An index in the range 0 ... (max number of entities in set - 1).
     */
    virtual IndexType entityIndex(const Entity<0>& e) const = 0;
    /** \brief Map entity of codimension 1 to index. */
    virtual IndexType entityIndex(const Entity<1>& e) const = 0;
    /** \brief Map entity of codimension 2 to index. */
    virtual IndexType entityIndex(const Entity<2>& e) const = 0;
    /** \brief Map entity of codimension 3 to index. */
    virtual IndexType entityIndex(const Entity<3>& e) const = 0;
};

/** \brief Wrapper of the index set specific to a Dune grid view class DuneGridView

 \internal The grid view class, rather than an index set class, is used as a
 template parameter because the latter doesn't provide information about the
 entity type.

 For consistency with IdSet it would be possible to take as parameters
 DuneGrid and DuneIndexSet instead.
 */
template<typename DuneGridView>
class ConcreteIndexSet: public IndexSet
{
public:
    typedef typename DuneGridView::IndexSet DuneIndexSet;

private:
    const DuneIndexSet* m_dune_index_set;

public:
    /**  Constructor */
    explicit ConcreteIndexSet(const DuneIndexSet* dune_index_set) :
        m_dune_index_set(dune_index_set) {
    }

    /** Read-only access to the underlying Dune index set */
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
    entityCodimNIndex(const Entity<codim>& e) const;
    template <int codim>
    typename boost::enable_if_c<codim <= DuneGridView::dimension, IndexType>::type
    entityCodimNIndex(const Entity<codim>& e) const;
};

} // namespace Bempp

#endif
