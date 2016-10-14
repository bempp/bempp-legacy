#ifndef bempp_p1_grid_hpp
#define bempp_p1_grid_hpp

#include "../../common/common.hpp"
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entityseed.hh>
#include <dune/grid/common/entitypointer.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/leveliterator.hh>
#include <boost/none.hpp> 


namespace BemppGrid {

    template <int mydim, int cdim, class> class P1GridGeometry;
    template <int> class P1EntityPointerImp;
    template <int> class P1EntitySeedImp;
    template <int codim, Dune::PartitionIteratorType, class> class P1LevelIteratorImp;

    struct P1GridFamily {};

    class P1Grid {

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
             typedef Dune::Geometry<2-cd, 3, P1Grid, P1GridGeometry> Geometry;
             typedef Dune::EntitySeed<P1Grid, P1EntitySeedImp<cd>>  EntitySeed;
             typedef boost::none_t LocalGeometry;
             typedef Dune::EntityPointer<P1Grid, P1EntityPointerImp<cd>> EntityPointer;
             typedef Dune::LevelIterator<cd, Dune::All_Partition, P1Grid, P1LevelIteratorImp> Iterator;

             template <Dune::PartitionIteratorType pitype>
             struct Partition {

                 typedef typename Dune::LevelIterator<cd, pitype, P1Grid, P1LevelIteratorImp> LevelIterator;
                 typedef typename Dune::LevelIterator<cd, pitype, P1Grid, P1LevelIteratorImp> LeafIterator;


             };             
             typedef Iterator LevelIterator;
             typedef Iterator LeafIterator;

            };

            typedef boost::none_t LeafIntersectionIterator;
            typedef boost::none_t LevelIntersectionIterator;
            typedef boost::none_t HierarchicIterator;


            // Codim< cd >::template Partition< pitype >::LevelIterator     lbegin (int level) const
    };

    

}

#include "p1_grid_imp.hpp"

#endif
