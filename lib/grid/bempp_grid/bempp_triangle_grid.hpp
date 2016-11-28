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
#include <dune/common/parallel/collectivecommunication.hh>
#include <boost/none.hpp> 


namespace BemppGrid {

    class TriangleGrid;

    class DataContainer;
    class LevelIndexSetImp;

    class P1GlobalIdSetImp;
    class P1LocalIdSetImp;


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
    template <class, Dune::PartitionIteratorType> class LevelGridViewTraits;
    template <class, Dune::PartitionIteratorType> class LeafGridViewTraits;

    struct TriangleGridFamily {

        typedef Dune::GridTraits<
            2, 3,
            TriangleGrid,
            Geometry,
            EntityImp,
            EntityPointerImp,
            LevelIteratorImp,
            P1LeafIntersectionImp,
            P1LevelIntersectionImp,
            P1LeafIntersectionIteratorImp,
            P1LevelIntersectionIteratorImp,
            P1HierarchicIteratorImp,
            LevelIteratorImp,
            LevelIndexSetImp,
            LevelIndexSetImp,
            P1GlobalIdSetImp,
            std::size_t,
            P1LocalIdSetImp,
            std::size_t,
            Dune::CollectiveCommunication<TriangleGrid>,
            LevelGridViewTraits,
            LeafGridViewTraits,
            EntitySeedImp>
                Traits;

        };

    class TriangleGrid 
        : public Dune::GridDefaultImplementation<2, 3, double, TriangleGridFamily> {

        private:

            friend class LevelIndexSetImp;
            friend class DataContainer;

            friend class EntityPointerImp<0, const TriangleGrid>;
            friend class EntityPointerImp<1, const TriangleGrid>;
            friend class EntityPointerImp<2, const TriangleGrid>;

        public: 

            typedef double ctype;
            typedef TriangleGridFamily GridFamily;

            enum {dimension = 2};
            enum {dimensionworld = 3};

            template <Dune::PartitionIteratorType pitype>
            struct Partition
            {
               typedef boost::none_t LevelGridView;
               typedef boost::none_t LeafGridView;
            };

            typedef typename Partition<Dune::All_Partition>::LevelGridView LevelGridView;
            typedef typename Partition<Dune::All_Partition>::LeafGridView LeafGridView;

            template <int cd>
            struct Codim
            {
             typedef Dune::Geometry<2-cd, 3, const TriangleGrid, Geometry> Geometry;
             typedef Dune::EntitySeed<const TriangleGrid, EntitySeedImp<cd, const TriangleGrid>>  EntitySeed;
             typedef boost::none_t LocalGeometry;
             typedef Dune::EntityPointer<const TriangleGrid, EntityPointerImp<cd, const TriangleGrid>> EntityPointer;
             typedef Dune::EntityIterator<cd, const TriangleGrid, LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>> Iterator;
             typedef Dune::Entity<cd, 2, const TriangleGrid, EntityImp> Entity;

             template <Dune::PartitionIteratorType pitype>
             struct Partition {

                 typedef typename Dune::EntityIterator<cd, const TriangleGrid, LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>> LevelIterator;
                 typedef typename Dune::EntityIterator<cd, const TriangleGrid, LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>> LeafIterator;


             };             
             typedef Iterator LevelIterator;
             typedef Iterator LeafIterator;

            };

            typedef boost::none_t LeafIntersectionIterator;
            typedef boost::none_t LevelIntersectionIterator;
            typedef boost::none_t HierarchicIterator;

            TriangleGrid(const shared_ptr<DataContainer>& data); 

            template <int cd>
            typename Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator lbegin(int level) const;

            template <int cd>
            typename Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator lend(int level) const;


        private:

            template <int cd>
            static GridFamily::Traits::LevelIndexSet::IndexType entityIndex(
                    const typename GridFamily::Traits::Codim<cd>::Entity& entity);

            template <int cd>
            static int entityLevel(
                    const typename GridFamily::Traits::Codim<cd>::Entity& entity);

            template <int cd>
            static const EntityImp<cd, 2, const TriangleGrid>& getEntityImp(
                    const typename GridFamily::Traits::Codim<cd>::Entity& entity);

            shared_ptr<DataContainer> m_data;
    };

    

}

#include "bempp_triangle_grid_imp.hpp"

#endif
