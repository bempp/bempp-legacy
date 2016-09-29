#ifndef bempp_p1_grid_hpp
#define bempp_p1_grid_hpp

#include "../../common/common.hpp"
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/entityseed.hh>
#include <boost/none.hpp>


namespace BemppGrid {

    template <int mydim, int cdim, class> class P1Geometry;
    template <int> class P1EntitySeedImp;

    struct P1GridFamily {};

    class P1Grid {

        public: 

            typedef double ctype;
            typedef P1GridFamily GridFamily;

            enum {dimension = 2};
            enum {dimensionworld = 3};

            template <int cd>
            struct Codim
            {
             typedef P1Geometry<2-cd, 3, P1Grid> Geometry;
             typedef Dune::EntitySeed<P1Grid, P1EntitySeedImp<cd>>  EntitySeed;
             typedef boost::none_t LocalGeometry;
             typedef boost::none_t EntityPointer;
            };

            typedef boost::none_t LeafIntersectionIterator;
            typedef boost::none_t LevelIntersectionIterator;
            typedef boost::none_t HierarchicIterator;

    
    };

    

}

#include "p1_impl/p1_entity_seed.hpp"
#include "p1_impl/p1_geometry.hpp"
#include "p1_impl/p1_entity.hpp"
#include "p1_impl/p1_entity_iterator.hpp"
#include "p1_impl/p1_data_container.hpp"

#endif
