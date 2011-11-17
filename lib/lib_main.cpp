/** \name Explicit instantiations of library classes @{ */

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"

namespace Bempp
{

// Grid
template class ConcreteGrid<DefaultDuneGrid> ;

// Grid views
template class ConcreteGridView<DefaultDuneGrid::LeafGridView> ;
template class ConcreteGridView<DefaultDuneGrid::LevelGridView> ;

/** \bug Why doesn't this work? */
// template std::auto_ptr<EntityIterator<0> >
// ConcreteGridView<DefaultDuneGrid::LeafGridView>::entityIterator<0>() const;
// template std::auto_ptr<EntityIterator<1> >
// ConcreteGridView<DefaultDuneGrid::LeafGridView>::entityIterator<1>() const;
// template std::auto_ptr<EntityIterator<2> >
// ConcreteGridView<DefaultDuneGrid::LeafGridView>::entityIterator<2>() const;
// template std::auto_ptr<EntityIterator<0> >
// ConcreteGridView<DefaultDuneGrid::LevelGridView>::entityIterator<0>() const;
// template std::auto_ptr<EntityIterator<1> >
// ConcreteGridView<DefaultDuneGrid::LevelGridView>::entityIterator<1>() const;
// template std::auto_ptr<EntityIterator<2> >
// ConcreteGridView<DefaultDuneGrid::LevelGridView>::entityIterator<2>() const;

// Id and index sets
template class ConcreteIdSet<DefaultDuneGrid, DefaultDuneGrid::GlobalIdSet> ;
template class ConcreteIndexSet<DefaultDuneGrid::LeafGridView> ;
template class ConcreteIndexSet<DefaultDuneGrid::LevelGridView> ;

// Entities
template class ConcreteEntity<0, DefaultDuneGrid::Codim<0>::Entity> ;
template class ConcreteEntity<1, DefaultDuneGrid::Codim<1>::Entity> ;
template class ConcreteEntity<2, DefaultDuneGrid::Codim<2>::Entity> ;

/** \bug Why doesn't this work? */
// template std::auto_ptr<EntityIterator<1> >
// ConcreteEntity<0, DefaultDuneGrid::Codim<0>::Entity>::subEntityIterator<1>() const;
// template std::auto_ptr<EntityIterator<2> >
// ConcreteEntity<0, DefaultDuneGrid::Codim<0>::Entity>::subEntityIterator<2>() const;

// Entity pointers
template class ConcreteEntityPointer<DefaultDuneGrid::Codim<0>::EntityPointer> ;
template class ConcreteEntityPointer<DefaultDuneGrid::Codim<1>::EntityPointer> ;
template class ConcreteEntityPointer<DefaultDuneGrid::Codim<2>::EntityPointer> ;

// Iterators
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<0>::LeafIterator> ;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<1>::LeafIterator> ;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<2>::LeafIterator> ;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<0>::LevelIterator> ;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<1>::LevelIterator> ;
template class ConcreteRangeEntityIterator<
DefaultDuneGrid::Codim<2>::LevelIterator> ;
template class ConcreteRangeEntityIterator<DefaultDuneGrid::HierarchicIterator> ;
template class ConcreteSubentityIterator<DefaultDuneGrid::Codim<0>::Entity, 1> ;
template class ConcreteSubentityIterator<DefaultDuneGrid::Codim<0>::Entity, 2> ;

// Geometries
template class ConcreteGeometry<DefaultDuneGrid::Codim<0>::Entity::Geometry> ;
template class ConcreteGeometry<DefaultDuneGrid::Codim<1>::Entity::Geometry> ;
template class ConcreteGeometry<DefaultDuneGrid::Codim<2>::Entity::Geometry> ;
}

/** @} */
