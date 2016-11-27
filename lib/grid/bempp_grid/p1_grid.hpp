#ifndef bempp_p1_grid_hpp
#define bempp_p1_grid_hpp

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

    class P1Grid;

    class P1DataContainer;
    class P1LevelIndexSetImp;

    class P1GlobalIdSetImp;
    class P1LocalIdSetImp;


    template <int mydim, int cdim, class> class Geometry;
    template <int codim, int dim, class> class P1EntityImp;
    template <int, class> class P1EntityPointerImp;
    template <int, class> class P1EntitySeedImp;
    template <int codim, Dune::PartitionIteratorType, class> class P1LevelIteratorImp;
    template <class> class P1LeafIntersectionImp;
    template <class> class P1LevelIntersectionImp;
    template <class> class P1HierarchicIteratorImp;
    template <class> class P1LeafIntersectionIteratorImp;
    template <class> class P1LevelIntersectionIteratorImp;
    template <class, Dune::PartitionIteratorType> class LevelGridViewTraits;
    template <class, Dune::PartitionIteratorType> class LeafGridViewTraits;

    struct P1GridFamily {

        typedef Dune::GridTraits<
            2, 3,
            P1Grid,
            Geometry,
            P1EntityImp,
            P1EntityPointerImp,
            P1LevelIteratorImp,
            P1LeafIntersectionImp,
            P1LevelIntersectionImp,
            P1LeafIntersectionIteratorImp,
            P1LevelIntersectionIteratorImp,
            P1HierarchicIteratorImp,
            P1LevelIteratorImp,
            P1LevelIndexSetImp,
            P1LevelIndexSetImp,
            P1GlobalIdSetImp,
            std::size_t,
            P1LocalIdSetImp,
            std::size_t,
            Dune::CollectiveCommunication<P1Grid>,
            LevelGridViewTraits,
            LeafGridViewTraits,
            P1EntitySeedImp>
                Traits;

        };

    class P1Grid 
        : public Dune::GridDefaultImplementation<2, 3, double, P1GridFamily> {

            friend class P1LevelIndexSetImp;

        public: 

            typedef double ctype;
            typedef P1GridFamily GridFamily;

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
             typedef Dune::Geometry<2-cd, 3, const P1Grid, Geometry> Geometry;
             typedef Dune::EntitySeed<const P1Grid, P1EntitySeedImp<cd, const P1Grid>>  EntitySeed;
             typedef boost::none_t LocalGeometry;
             typedef Dune::EntityPointer<const P1Grid, P1EntityPointerImp<cd, const P1Grid>> EntityPointer;
             typedef Dune::EntityIterator<cd, const P1Grid, P1LevelIteratorImp<cd, Dune::All_Partition, const P1Grid>> Iterator;
             typedef Dune::Entity<cd, 2, const P1Grid, P1EntityImp> Entity;

             template <Dune::PartitionIteratorType pitype>
             struct Partition {

                 typedef typename Dune::EntityIterator<cd, const P1Grid, P1LevelIteratorImp<cd, Dune::All_Partition, const P1Grid>> LevelIterator;
                 typedef typename Dune::EntityIterator<cd, const P1Grid, P1LevelIteratorImp<cd, Dune::All_Partition, const P1Grid>> LeafIterator;


             };             
             typedef Iterator LevelIterator;
             typedef Iterator LeafIterator;

            };

            typedef boost::none_t LeafIntersectionIterator;
            typedef boost::none_t LevelIntersectionIterator;
            typedef boost::none_t HierarchicIterator;

            P1Grid(const shared_ptr<P1DataContainer>& data); 

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

            shared_ptr<P1DataContainer> m_data;
    };

    

}

#include "p1_grid_imp.hpp"

#endif
