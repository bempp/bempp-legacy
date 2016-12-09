#ifndef bempp_grid_triangle_imp_entity_hpp
#define bempp_grid_triangle_imp_entity_hpp

#include "../../../common/common.hpp"
#include "../bempp_grid_types.hpp"
#include "bempp_grid_geometry.hpp"
#include "bempp_grid_entity_seed.hpp"
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/entityseed.hh>
#include <dune/grid/common/entitypointer.hh>

namespace BemppGrid {

class TriangleGrid;
class DataContainer;
class LevelIndexSetImp;

template <int cd, int dim, class GridImp> class EntityImp {};

template <int codim, class> class EntityPointerImp;

// Elements
template <>
class EntityImp<0, 2, const TriangleGrid>
    : public Dune::EntityDefaultImplementation<0, 2, const TriangleGrid,
                                               EntityImp> {

  friend class EntitySeedImp<0, const TriangleGrid>;
  friend class EntityPointerImp<0, const TriangleGrid>;
  friend class TriangleGrid;

public:
  typedef
      typename TriangleGrid::GridFamily::Traits::Codim<0>::Geometry Geometry;
  typedef typename TriangleGrid::GridFamily::Traits::Codim<0>::Entity Entity;
  typedef typename TriangleGrid::GridFamily::Traits::Codim<0>::EntitySeed
      EntitySeed;

  template <int codim>
  using EntityPointer =
      Dune::EntityPointer<const TriangleGrid,
                          EntityPointerImp<codim, const TriangleGrid>>;

  EntityImp(const shared_ptr<DataContainer> &data, int level,
            unsigned int index);

  int level() const;
  Geometry geometry() const;
  Dune::GeometryType type() const;
  Dune::PartitionType partitionType() const;
  bool equals(const EntityImp<0, 2, const TriangleGrid> &other) const;
  EntitySeed seed() const;

  template <int cd> EntityPointer<cd> subEntity(int i) const;

  template <int cc> int count() const;

  bool hasFather() const;
  bool isLeaf() const;
  bool isRegular() const;
  bool isNew() const;

  EntityPointer<0> father() const;

  unsigned int id() const;

  unsigned int subEntities(unsigned int codim) const;

private:
  EntityImp() {}

  shared_ptr<DataContainer> m_data;
  int m_level;
  unsigned int m_index;
};

// Edges
template <>
class EntityImp<1, 2, const TriangleGrid>
    : public Dune::EntityDefaultImplementation<1, 2, const TriangleGrid,
                                               EntityImp> {

  friend class EntityPointerImp<1, const TriangleGrid>;
  friend class TriangleGrid;

public:
  typedef TriangleGrid::GridFamily::Traits::Codim<1>::EntitySeed EntitySeed;
  typedef
      typename TriangleGrid::GridFamily::Traits::Codim<1>::Geometry Geometry;
  typedef typename TriangleGrid::GridFamily::Traits::Codim<1>::Entity Entity;

  EntityImp(const shared_ptr<DataContainer> &data, int level,
            unsigned int index);

  int level() const;
  Geometry geometry() const;
  Dune::GeometryType type() const;
  Dune::PartitionType partitionType() const;
  bool equals(const EntityImp<1, 2, const TriangleGrid> &other) const;
  EntitySeed seed() const;

  unsigned int id() const;

private:
  EntityImp() {}

  shared_ptr<DataContainer> m_data;
  int m_level;
  unsigned int m_index;
};

// Vertices
template <>
class EntityImp<2, 2, const TriangleGrid>
    : public Dune::EntityDefaultImplementation<2, 2, const TriangleGrid,
                                               EntityImp> {

  friend class EntityPointerImp<2, const TriangleGrid>;
  friend class TriangleGrid;

public:
  typedef TriangleGrid::GridFamily::Traits::Codim<2>::EntitySeed EntitySeed;
  typedef
      typename TriangleGrid::GridFamily::Traits::Codim<2>::Geometry Geometry;
  typedef typename TriangleGrid::GridFamily::Traits::Codim<2>::Entity Entity;

  EntityImp(const shared_ptr<DataContainer> &data, int level,
            unsigned int index);

  int level() const;
  Geometry geometry() const;
  Dune::GeometryType type() const;
  Dune::PartitionType partitionType() const;
  bool equals(const EntityImp<2, 2, const TriangleGrid> &other) const;
  EntitySeed seed() const;

  unsigned int id() const;

private:
  EntityImp() {}

  shared_ptr<DataContainer> m_data;
  int m_level;
  unsigned int m_index;
};
}

#include "bempp_grid_entity_imp.hpp"

#endif
