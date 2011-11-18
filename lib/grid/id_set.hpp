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

#ifndef bempp_id_set_hpp
#define bempp_id_set_hpp

#include "id_set_decl.hpp"
#include "entity.hpp"

#include <stdexcept>

namespace Bempp
{

template<typename DuneGrid, typename DuneIdSet>
template<int codim>
inline typename boost::disable_if_c<(codim <= DuneGrid::dimension), IdSet::IdType>::type
ConcreteIdSet<DuneGrid, DuneIdSet>::entityCodimNId(const Entity<codim>& e) const
{
    throw std::logic_error("IdSet::entityId(): invalid entity codimension");
}

template<typename DuneGrid, typename DuneIdSet>
template<int codim>
inline typename boost::enable_if_c<(codim <= DuneGrid::dimension), IdSet::IdType>::type
ConcreteIdSet<DuneGrid, DuneIdSet>::entityCodimNId(const Entity<codim>& e) const
{
    typedef typename DuneGrid::template Codim<codim>::Entity DuneEntity;
    typedef ConcreteEntity<codim, DuneEntity> ConcEntity;
    const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
    return m_dune_id_set->id(ce.duneEntity());
}

} // namespace Bempp

#endif // ID_SET_HPP
