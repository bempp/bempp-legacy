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

#ifndef bempp_concrete_element_mapper_hpp
#define bempp_concrete_element_mapper_hpp

#include "mapper.hpp"

#include <dune/grid/common/mcmgmapper.hh>

namespace Bempp
{

/** \brief An injective mapping from the full set of codimension-0 entities
  ("elements") of a grid to the integers 0 ... (number of entities - 1). */
template <typename DuneGridView>
class ConcreteElementMapper : public Mapper
{
private:
    Dune::MultipleCodimMultipleGeomTypeMapper<DuneGridView,
    Dune::MCMGElementLayout > m_dune_mapper;

public:
    ConcreteElementMapper(const DuneGridView& dune_gv) :
        m_dune_mapper(dune_gv,
                      Dune::MCMGElementLayout<DuneGridView::dimension>())
    {}

    virtual int size() const {
        return m_dune_mapper.size();
    }

    virtual int entityIndex(const Entity<0>& e) const {
        typedef typename DuneGridView::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_mapper.map(ce.duneEntity());
    }

    virtual int entityIndex(const Entity<1>& e) const {
        throw std::logic_error("ConcreteElementMapper::entityIndex(): "
                               "entities of codimension 1 do not belong to the "
                               "managed set.");
    }

    virtual int entityIndex(const Entity<2>& e) const {
        throw std::logic_error("ConcreteElementMapper::entityIndex(): "
                               "entities of codimension 2 do not belong to the "
                               "managed set.");
    }

    virtual int entityIndex(const Entity<3>& e) const {
        throw std::logic_error("ConcreteElementMapper::entityIndex(): "
                               "entities of codimension 2 do not belong to the "
                                "managed set.");
    }

    virtual int subEntityIndex(const Entity<0>& e, int i,
                               unsigned int codimSub) const {
        typedef typename DuneGridView::template Codim<0>::Entity DuneEntity;
        typedef ConcreteEntity<0, DuneEntity> ConcEntity;
        const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
        return m_dune_mapper.map(ce.duneEntity(), i, codimSub);
    }
};

} // namespace Bempp

#endif
