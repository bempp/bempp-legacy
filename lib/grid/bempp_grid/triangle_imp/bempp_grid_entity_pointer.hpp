#ifndef bempp_grid_triangle_imp_entity_pointer_hpp
#define bempp_grid_triangle_imp_entity_pointer_hpp

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/gridenums.hh>
#include <memory>

namespace BemppGrid {

class TriangleGrid;

template <int, int, class> class EntityImp;
template <int, Dune::PartitionIteratorType, class> class LevelIteratorImp;

template <int codim, class> class EntityPointerImp {

  template <int cd, Dune::PartitionIteratorType pitype, class GridType>
  friend class LevelIteratorImp;

public:
  enum { codimension = codim };
  typedef
      typename TriangleGrid::GridFamily::Traits::Codim<codim>::Entity Entity;

  EntityPointerImp(const EntityImp<codim, 2, const TriangleGrid> &entity)
      : m_entity(new Entity(entity)) {}

  EntityPointerImp(const EntityPointerImp<codim, const TriangleGrid> &other)
      : m_entity(
            new Entity(TriangleGrid::getEntityImp<codim>(*other.m_entity))) {}

  EntityPointerImp<codim, const TriangleGrid> &
  operator=(const EntityPointerImp<codim, const TriangleGrid> &other) {
    setEntity(TriangleGrid::getEntityImp<codim>(*other.m_entity));
    return *this;
  }

  virtual ~EntityPointerImp() {}

  Entity &dereference() const { return *m_entity; }

  int level() const { return m_entity->level(); }

  bool equals(const EntityPointerImp<codim, const TriangleGrid> &other) const {

    return (other.level() == level()) &&
           (TriangleGrid::entityIndex<codim>(*other.m_entity) ==
            TriangleGrid::entityIndex<codim>(*m_entity));
  }

private:
  void setEntity(const EntityImp<codim, 2, const TriangleGrid> &entity) const {

    // Need pointer as Entity has no assignment operator
    m_entity.reset(new Entity(entity));
  }

  mutable std::unique_ptr<Entity> m_entity;
};
}

#endif
