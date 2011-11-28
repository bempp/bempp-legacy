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

/** \file Explicit instantiations of library classes.
    Needed for Python and Matlab support.
*/

#include "entity.hpp"
#include "entity_iterator.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "grid_view.hpp"

namespace Bempp
{

// Grid
template class ConcreteGrid<DefaultDuneGrid>;

// Grid views
template class ConcreteGridView<DefaultDuneGrid::LeafGridView>;
template class ConcreteGridView<DefaultDuneGrid::LevelGridView>;

template std::auto_ptr<EntityIterator<0> > GridView::entityIterator<0>() const;
template std::auto_ptr<EntityIterator<1> > GridView::entityIterator<1>() const;
template std::auto_ptr<EntityIterator<2> > GridView::entityIterator<2>() const;
template std::auto_ptr<EntityIterator<3> > GridView::entityIterator<3>() const;

// Id and index sets
template class ConcreteIdSet<DefaultDuneGrid, DefaultDuneGrid::GlobalIdSet>;
template class ConcreteIndexSet<DefaultDuneGrid::LeafGridView>;
template class ConcreteIndexSet<DefaultDuneGrid::LevelGridView>;

// Entities
template class ConcreteEntity<0, DefaultDuneGrid::Codim<0>::Entity>;
template class ConcreteEntity<1, DefaultDuneGrid::Codim<1>::Entity>;
template class ConcreteEntity<2, DefaultDuneGrid::Codim<2>::Entity>;

template std::auto_ptr<EntityIterator<1> > Entity<0>::subEntityIterator<1>() const;
template std::auto_ptr<EntityIterator<2> > Entity<0>::subEntityIterator<2>() const;
template std::auto_ptr<EntityIterator<3> > Entity<0>::subEntityIterator<3>() const;

// Entity pointers
template class ConcreteEntityPointer<DefaultDuneGrid::Codim<0>::EntityPointer>;
template class ConcreteEntityPointer<DefaultDuneGrid::Codim<1>::EntityPointer>;
template class ConcreteEntityPointer<DefaultDuneGrid::Codim<2>::EntityPointer>;

// Iterators
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<0>::LeafIterator>;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<1>::LeafIterator>;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<2>::LeafIterator>;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<0>::LevelIterator>;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<1>::LevelIterator>;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<2>::LevelIterator>;
template class ConcreteRangeEntityIterator<DefaultDuneGrid::HierarchicIterator>;
template class ConcreteSubentityIterator<DefaultDuneGrid::Codim<0>::Entity, 1>;
template class ConcreteSubentityIterator<DefaultDuneGrid::Codim<0>::Entity, 2>;

// Geometries
template class ConcreteGeometry<DefaultDuneGrid::Codim<0>::Entity::Geometry>;
template class ConcreteGeometry<DefaultDuneGrid::Codim<1>::Entity::Geometry>;
template class ConcreteGeometry<DefaultDuneGrid::Codim<2>::Entity::Geometry>;

// VtkWriters

template class ConcreteVtkWriter<DefaultDuneGrid::LeafGridView>;
template class ConcreteVtkWriter<DefaultDuneGrid::LevelGridView>;

}
