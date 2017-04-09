#ifndef bempp_triangle_grid_hpp
#define bempp_triangle_grid_hpp

#include "../../common/common.hpp"
#include "bempp_grid_types.hpp"
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entityseed.hh>
#include <dune/grid/common/entitypointer.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/entityiterator.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/defaultgridview.hh>
#include <dune/geometry/type.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <boost/none.hpp>

namespace BemppGrid {

class TriangleGrid;

class DataContainer;
class LevelIndexSetImp;

class IdSetImp;

template <int mydim, int cdim, class> class Geometry;
template <int codim, int dim, class> class EntityImp;
template <int, class> class EntityPointerImp;
template <int, class> class EntitySeedImp;
template <int codim, Dune::PartitionIteratorType, class> class LevelIteratorImp;
template <class> class P1LeafIntersectionImp;
template <class> class P1LevelIntersectionImp;
template <class> class P1HierarchicIteratorImp;
template <class> class P1LeafIntersectionIteratorImp;
template <class> class P1LevelIntersectionIteratorImp;

struct TriangleGridFamily {

#if DUNE_VERSION <= 241

  typedef Dune::GridTraits<
      2, 3, TriangleGrid, Geometry, EntityImp, EntityPointerImp,
      LevelIteratorImp, P1LeafIntersectionImp, P1LevelIntersectionImp,
      P1LeafIntersectionIteratorImp, P1LevelIntersectionIteratorImp,
      P1HierarchicIteratorImp, LevelIteratorImp, LevelIndexSetImp,
      LevelIndexSetImp, IdSetImp, unsigned int, IdSetImp, unsigned int,
      Dune::CollectiveCommunication<TriangleGrid>,
      Dune::DefaultLevelGridViewTraits, Dune::DefaultLeafGridViewTraits,
      EntitySeedImp> Traits;

#else 

      typedef Dune::GridTraits<
          2, 3, TriangleGrid, Geometry, EntityImp,
          LevelIteratorImp, P1LeafIntersectionImp, P1LevelIntersectionImp,
          P1LeafIntersectionIteratorImp, P1LevelIntersectionIteratorImp,
          P1HierarchicIteratorImp, LevelIteratorImp, LevelIndexSetImp,
          LevelIndexSetImp, IdSetImp, unsigned int, IdSetImp, unsigned int,
          Dune::CollectiveCommunication<TriangleGrid>,
          Dune::DefaultLevelGridViewTraits, Dune::DefaultLeafGridViewTraits,
          EntitySeedImp> Traits;
#endif

};

class TriangleGrid
    : public Dune::GridDefaultImplementation<2, 3, double, TriangleGridFamily> {

private:
  friend class LevelIndexSetImp;
  friend class IdSetImp;
  friend class DataContainer;

  friend class EntityPointerImp<0, const TriangleGrid>;
  friend class EntityPointerImp<1, const TriangleGrid>;
  friend class EntityPointerImp<2, const TriangleGrid>;

public:
  typedef double ctype;
  typedef TriangleGridFamily GridFamily;

  enum { dimension = 2 };
  enum { dimensionworld = 3 };

#if DUNE_VERSION <= 241

  template <Dune::PartitionIteratorType pitype> struct Partition {
    typedef typename TriangleGridFamily::Traits::template Partition<
        pitype>::LevelGridView LevelGridView;
    typedef typename TriangleGridFamily::Traits::template Partition<
        pitype>::LeafGridView LeafGridView;
  };

  typedef typename Partition<Dune::All_Partition>::LevelGridView LevelGridView;
  typedef typename Partition<Dune::All_Partition>::LeafGridView LeafGridView;

#else

  typedef typename TriangleGridFamily::Traits::LevelGridView LevelGridView;
  typedef typename TriangleGridFamily::Traits::LeafGridView LeafGridView;

#endif

  template <int cd> struct Codim {
    typedef Dune::Geometry<2 - cd, 3, const TriangleGrid, Geometry> Geometry;
    typedef Dune::EntitySeed<const TriangleGrid,
                             EntitySeedImp<cd, const TriangleGrid>> EntitySeed;
    typedef boost::none_t LocalGeometry;
    typedef Dune::EntityPointer<const TriangleGrid,
                                EntityPointerImp<cd, const TriangleGrid>>
        EntityPointer;
    typedef Dune::EntityIterator<
        cd, const TriangleGrid,
        LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>> Iterator;
    typedef Dune::Entity<cd, 2, const TriangleGrid, EntityImp> Entity;

    template <Dune::PartitionIteratorType pitype> struct Partition {

      typedef typename Dune::EntityIterator<
          cd, const TriangleGrid,
          LevelIteratorImp<cd, pitype, const TriangleGrid>> LevelIterator;
      typedef typename Dune::EntityIterator<
          cd, const TriangleGrid,
          LevelIteratorImp<cd, pitype, const TriangleGrid>> LeafIterator;
    };
    typedef Iterator LevelIterator;
    typedef Iterator LeafIterator;
  };

  typedef boost::none_t LeafIntersectionIterator;
  typedef boost::none_t LevelIntersectionIterator;
  typedef boost::none_t HierarchicIterator;

  TriangleGrid(const shared_ptr<DataContainer> &data);

  template <int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::LevelIterator
  lbegin(int level) const;

  template <int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::LevelIterator
  lend(int level) const;

  template <int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::LevelIterator
  leafbegin() const;

  template <int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::LevelIterator leafend() const;

  const typename GridFamily::Traits::LevelIndexSet &
  levelIndexSet(int level) const;
  const typename GridFamily::Traits::LeafIndexSet &leafIndexSet() const;

  const typename GridFamily::Traits::GlobalIdSet &globalIdSet() const;
  const typename GridFamily::Traits::LocalIdSet &localIdSet() const;

  int maxLevel() const;
  int size(int level, int codim) const;
  int size(int codim) const;
  int size(int level, Dune::GeometryType type) const;
  int size(Dune::GeometryType type) const;

private:
  const DataContainer &data() const;

  template <int cd>
  static GridFamily::Traits::LevelIndexSet::IndexType
  entityIndex(const typename GridFamily::Traits::Codim<cd>::Entity &entity);

  template <int cd>
  static int
  entityLevel(const typename GridFamily::Traits::Codim<cd>::Entity &entity);

  template <int cd>
  static const EntityImp<cd, 2, const TriangleGrid> &
  getEntityImp(const typename GridFamily::Traits::Codim<cd>::Entity &entity);

  shared_ptr<DataContainer> m_data;

  std::vector<shared_ptr<typename GridFamily::Traits::LevelIndexSet>>
      m_levelIndexSet;

  const typename GridFamily::Traits::GlobalIdSet m_globalIdSet;
  const typename GridFamily::Traits::LocalIdSet m_localIdSet;
};
}

#include "bempp_triangle_grid_imp.hpp"

#endif
