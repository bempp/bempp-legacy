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

#ifndef bempp_lib_grid_grid_hpp
#define bempp_lib_grid_grid_hpp

#include "grid_decl.hpp"
#include "entity.hpp"
#include "grid_view.hpp"
#include "id_set.hpp"

#include <stack> // fix a bug in foamgrid -- this header is not included where it should be
#include <dune/foamgrid/foamgrid.hh>

namespace Bempp
{

template<typename DuneGrid>
std::auto_ptr<GridView> ConcreteGrid<DuneGrid>::levelView(int level) const
{
    return std::auto_ptr<GridView>(new ConcreteGridView<typename DuneGrid::LevelGridView>(
                                       m_dune_grid->levelView(level)));
}

template<typename DuneGrid>
std::auto_ptr<GridView> ConcreteGrid<DuneGrid>::leafView() const
{
    return std::auto_ptr<GridView>(new ConcreteGridView<typename DuneGrid::LeafGridView>(
                                       m_dune_grid->leafView()));
}

template<typename DuneGrid>
inline bool ConcreteGrid<DuneGrid>::mark(int refCount, const Entity<0>& e)
{
    /// FIXME: should we catch std::bad_cast or leave it to the user?
    typedef typename DuneGrid::template Codim<0>::Entity DuneEntity;
    typedef ConcreteEntity<0, DuneEntity> ConcEntity;
    const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
    return m_dune_grid->mark(refCount, ce.duneEntity());
}

template<typename DuneGrid>
inline int ConcreteGrid<DuneGrid>::getMark(const Entity<0>& e) const
{
    /// FIXME: should we catch std::bad_cast or leave it to the user?
    typedef typename DuneGrid::template Codim<0>::Entity DuneEntity;
    typedef ConcreteEntity<0, DuneEntity> ConcEntity;
    const ConcEntity& ce = dynamic_cast<const ConcEntity&>(e);
    return m_dune_grid->getMark(ce.duneEntity());
}

// Default grid typedefs
typedef Dune::FoamGrid<3> DefaultDuneGrid; // 3 -> dimWorld
typedef ConcreteGrid<DefaultDuneGrid> DefaultGrid;

} // namespace Bempp

#endif // GRID_HPP
