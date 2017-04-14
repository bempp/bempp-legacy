#ifndef bempp_grid_triangle_imp_entity_imp_hpp
#define bempp_grid_triangle_imp_entity_imp_hpp

#include "bempp_grid_entity.hpp"
#include "bempp_grid_data_container.hpp"
#include <dune/geometry/type.hh>

namespace BemppGrid {

inline EntityImp<0, 2, const TriangleGrid>::EntityImp(
    const shared_ptr<DataContainer> &data, int level, unsigned int index)
    : m_data(data), m_level(level), m_index(index) {}

inline EntityImp<0, 2, const TriangleGrid>::Geometry
EntityImp<0, 2, const TriangleGrid>::geometry() const {

  return m_data->geometry<0>(Entity<0>(*this));
}

inline Dune::PartitionType
EntityImp<0, 2, const TriangleGrid>::partitionType() const {

  return Dune::InteriorEntity;
}

inline int EntityImp<0, 2, const TriangleGrid>::level() const {
  return m_level;
}

inline Dune::GeometryType EntityImp<0, 2, const TriangleGrid>::type() const {

  Dune::GeometryType geometryType;
  geometryType.makeTriangle();
  return geometryType;
}

inline bool EntityImp<0, 2, const TriangleGrid>::equals(
    const EntityImp<0, 2, const TriangleGrid> &other) const {

  return m_level == other.m_level && m_index == other.m_index;
}

inline EntityImp<0, 2, const TriangleGrid>::EntitySeed
EntityImp<0, 2, const TriangleGrid>::seed() const {

  return EntityImp<0, 2, const TriangleGrid>::EntitySeed(
      EntitySeedImp<0, const TriangleGrid>(m_level, m_index));
}

template <>
inline EntityImp<0, 2, const TriangleGrid>::Entity<1>
EntityImp<0, 2, const TriangleGrid>::subEntity<1>(int i) const {

  return Entity<1>(EntityImp<1, 2, const TriangleGrid>(
      m_data, m_level, m_data->element2Edges(m_level, m_index)[i]));
}

template <>
inline EntityImp<0, 2, const TriangleGrid>::Entity<2>
EntityImp<0, 2, const TriangleGrid>::subEntity<2>(int i) const {

  return Entity<2>(EntityImp<2, 2, const TriangleGrid>(
      m_data, m_level, m_data->elements(m_level)[m_index][i]));
}

template <>
inline EntityImp<0, 2, const TriangleGrid>::Entity<0>
EntityImp<0, 2, const TriangleGrid>::subEntity<0>(int i) const {

  return Entity<0>(*this);
}

inline EntityImp<0, 2, const TriangleGrid>::Entity<0>
EntityImp<0, 2, const TriangleGrid>::father() const {

  int fatherIndex = m_data->getElementFatherIndex(m_level, m_index);
  return Entity<0>(
      EntityImp<0, 2, const TriangleGrid>(m_data, m_level, fatherIndex));
}

template <> inline int EntityImp<0, 2, const TriangleGrid>::count<0>() const {

  return 1;
}

template <> inline int EntityImp<0, 2, const TriangleGrid>::count<1>() const {

  return 3;
}

template <> inline int EntityImp<0, 2, const TriangleGrid>::count<2>() const {

  return 3;
}

inline bool EntityImp<0, 2, const TriangleGrid>::hasFather() const {

  return (m_level > 0);
}

inline bool EntityImp<0, 2, const TriangleGrid>::isLeaf() const {

  return (m_level == m_data->levels() - 1);
}

inline bool EntityImp<0, 2, const TriangleGrid>::isRegular() const {

  return true;
}

inline bool EntityImp<0, 2, const TriangleGrid>::isNew() const {

  if (m_level == 0)
    return false;

  int fatherIndex = m_data->getElementFatherIndex(m_level, m_index);
  return m_data->getElementSons(m_level - 1, fatherIndex).size() > 0;
}

inline unsigned int
EntityImp<0, 2, const TriangleGrid>::subEntities(unsigned int codim) const {

  if (codim == 0)
    return 1;

  if (codim == 1)
    return 3;

  if (codim == 2)
    return 3;

  throw std::runtime_error(
      "Entity::subEntities(): Error. Require 0 <= codim <= 2.");
}

inline EntityImp<1, 2, const TriangleGrid>::EntityImp(
    const shared_ptr<DataContainer> &data, int level, unsigned int index)
    : m_data(data), m_level(level), m_index(index) {}

inline EntityImp<1, 2, const TriangleGrid>::Geometry
EntityImp<1, 2, const TriangleGrid>::geometry() const {

  return m_data->geometry<1>(Entity<1>(*this));
}

inline int EntityImp<1, 2, const TriangleGrid>::level() const {
  return m_level;
}

inline Dune::GeometryType EntityImp<1, 2, const TriangleGrid>::type() const {

  Dune::GeometryType geometryType;
  geometryType.makeLine();
  return geometryType;
}

inline Dune::PartitionType
EntityImp<1, 2, const TriangleGrid>::partitionType() const {

  return Dune::InteriorEntity;
}

inline bool EntityImp<1, 2, const TriangleGrid>::equals(
    const EntityImp<1, 2, const TriangleGrid> &other) const {

  return m_level == other.m_level && m_index == other.m_index;
}

inline EntityImp<1, 2, const TriangleGrid>::EntitySeed
EntityImp<1, 2, const TriangleGrid>::seed() const {

  return EntityImp<1, 2, const TriangleGrid>::EntitySeed(
      EntitySeedImp<1, const TriangleGrid>(m_level, m_index));
}

template <>
inline EntityImp<1, 2, const TriangleGrid>::Entity<0>
EntityImp<1, 2, const TriangleGrid>::subEntity<0>(int i) const {

  throw std::runtime_error("Subentities of edges not defined for codim 0.");
}

template <>
inline EntityImp<1, 2, const TriangleGrid>::Entity<2>
EntityImp<1, 2, const TriangleGrid>::subEntity<2>(int i) const {

  return Entity<2>(EntityImp<2, 2, const TriangleGrid>(
      m_data, m_level, m_data->edges(m_level)[m_index][i]));
}



template <>
inline EntityImp<1, 2, const TriangleGrid>::Entity<1>
EntityImp<1, 2, const TriangleGrid>::subEntity<1>(int i) const {

  return Entity<1>(*this);
}


inline EntityImp<2, 2, const TriangleGrid>::EntityImp(
    const shared_ptr<DataContainer> &data, int level, unsigned int index)
    : m_data(data), m_level(level), m_index(index) {}

inline EntityImp<2, 2, const TriangleGrid>::EntitySeed
EntityImp<2, 2, const TriangleGrid>::seed() const {

  return EntityImp<2, 2, const TriangleGrid>::EntitySeed(
      EntitySeedImp<2, const TriangleGrid>(m_level, m_index));
}

inline EntityImp<2, 2, const TriangleGrid>::Geometry
EntityImp<2, 2, const TriangleGrid>::geometry() const {

  return m_data->geometry<2>(Entity<2>(*this));
}

inline int EntityImp<2, 2, const TriangleGrid>::level() const {
  return m_level;
}

inline Dune::GeometryType EntityImp<2, 2, const TriangleGrid>::type() const {

  Dune::GeometryType geometryType;
  geometryType.makeVertex();
  return geometryType;
}

inline Dune::PartitionType
EntityImp<2, 2, const TriangleGrid>::partitionType() const {

  return Dune::InteriorEntity;
}

inline bool EntityImp<2, 2, const TriangleGrid>::equals(
    const EntityImp<2, 2, const TriangleGrid> &other) const {

  return m_level == other.m_level && m_index == other.m_index;
}

template <>
inline EntityImp<2, 2, const TriangleGrid>::Entity<0>
EntityImp<2, 2, const TriangleGrid>::subEntity<0>(int i) const {

  throw std::runtime_error("Subentities of vertices not defined for codim 0.");
}

template <>
inline EntityImp<2, 2, const TriangleGrid>::Entity<1>
EntityImp<2, 2, const TriangleGrid>::subEntity<1>(int i) const {

  throw std::runtime_error("Subentities of vertices not defined for codim 1.");
}


template <>
inline EntityImp<2, 2, const TriangleGrid>::Entity<2>
EntityImp<2, 2, const TriangleGrid>::subEntity<2>(int i) const {

  return Entity<2>(*this);
}


inline unsigned int EntityImp<0, 2, const TriangleGrid>::id() const {

  return m_data->id<0>(Entity<0>(*this));
}

inline unsigned int EntityImp<1, 2, const TriangleGrid>::id() const {

  return m_data->id<1>(Entity<1>(*this));
}

inline unsigned int EntityImp<2, 2, const TriangleGrid>::id() const {

  return m_data->id<2>(Entity<2>(*this));
}
}

#endif
