#ifndef bempp_grid_triangle_imp_level_iterator_hpp
#define bempp_grid_triangle_imp_level_iterator_hpp

#include "bempp_grid_entity_pointer.hpp"
#include <dune/grid/common/gridenums.hh>

namespace BemppGrid {

class DataContainer;
class TriangleGrid;

template <int, int, class> class EntityImp;

template <int codim, Dune::PartitionIteratorType, class>
class LevelIteratorImp : public EntityPointerImp<codim, const TriangleGrid> {

public:
  typedef Dune::Entity<codim, 2, const TriangleGrid, EntityImp> Entity;
  typedef EntityPointerImp<codim, const TriangleGrid> Base;

  LevelIteratorImp(const shared_ptr<DataContainer> &data, int level,
                   unsigned int index)
      : m_data(data), m_level(level), m_index(index),
        Base(EntityPointerImp<codim, const TriangleGrid>(
            EntityImp<codim, 2, const TriangleGrid>(data, level, index))){};

  void increment() const {

    m_index++;
    static_cast<const Base *>(this)->setEntity(
        EntityImp<codim, 2, const TriangleGrid>(m_data, m_level, m_index));
  }

private:
  shared_ptr<DataContainer> m_data;
  int m_level;
  mutable unsigned int m_index;
};
}

#endif
