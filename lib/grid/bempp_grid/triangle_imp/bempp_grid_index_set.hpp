#ifndef bempp_grid_triangle_imp_index_set_hpp
#define bempp_grid_triangle_imp_index_set_hpp

#include "../../../common/common.hpp"
#include "bempp_grid_data_container.hpp"
#include <dune/grid/common/entity.hh>
#include <dune/geometry/type.hh>
#include <string>

namespace BemppGrid {

class TriangleGrid;
class DataContainer;
template <int, int, class> class EntityImp;

class LevelIndexSetImp
    : public Dune::IndexSetDefaultImplementation<const TriangleGrid,
                                                 LevelIndexSetImp> {

public:
  typedef unsigned int IndexType;

  LevelIndexSetImp(const shared_ptr<DataContainer> &data, unsigned int level)
      : m_data(data), m_level(level),
        m_types({{Dune::GeometryType(Dune::GeometryType::simplex, 2)},
                 {Dune::GeometryType(Dune::GeometryType::simplex, 1)},
                 {Dune::GeometryType(Dune::GeometryType::simplex, 0)}}) {}

  template <int cd>
  IndexType
  index(const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity
            &entity) const {

    return TriangleGrid::entityIndex<cd>(entity);
  }

  template <int cd>
  IndexType
  subIndex(const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity
               &entity,
           int i, unsigned int codim) const {

    if (codim < cd)
      throw std::runtime_error("LevelIndexSetImp::subIndex(): Error: cd = " +
                               std::to_string(cd) + ", codim = " +
                               std::to_string(codim));

    const auto& entityImpl = TriangleGrid::getEntityImp<cd>(entity);

    if (codim == 0)
      return index<0>(entityImpl.template subEntity<0>(i));
    if (codim == 1)
      return index<1>(entityImpl.template subEntity<1>(i));
    if (codim == 2)
      return index<2>(entityImpl.template subEntity<2>(i));
    throw std::runtime_error(
        "LevelIndexSetImp::subIndex(): Require 0 <= codim <= 2");
  }

  const std::vector<Dune::GeometryType> &geomTypes(int codim) const {
    return m_types[codim];
  }

  IndexType size(Dune::GeometryType type) const {

    if (type.isVertex())
      return m_data->numberOfEntities<2>(m_level);
    if (type.isLine())
      return m_data->numberOfEntities<1>(m_level);
    if (type.isTriangle())
      return m_data->numberOfEntities<0>(m_level);

    throw std::runtime_error(
        "LevelIndexSetImp::size(): Unknownn Geometry type");
  }

  IndexType size(int codim) const {

    if (codim == 2)
      return m_data->numberOfEntities<2>(m_level);
    if (codim == 1)
      return m_data->numberOfEntities<1>(m_level);
    if (codim == 0)
      return m_data->numberOfEntities<0>(m_level);

    throw std::runtime_error(
        "LevelIndexSetImp::size(): Unknownn Geometry type");
  }

  template <class EntityType> bool contains(const EntityType &e) const {

    return e.level() == m_level;
  }

private:
  shared_ptr<DataContainer> m_data;
  std::vector<std::vector<Dune::GeometryType>> m_types;
  unsigned int m_level;
};
}

#endif
