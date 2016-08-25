#ifndef bempp_p1_grid_hpp
#define bempp_p1_grid_hpp

#include "../../common/common.hpp"
#include "p1_impl/p1_geometry.hpp"
#include "p1_impl/p1_entity.hpp"
#include "p1_impl/p1_entity_iterator.hpp"
#include "p1_impl/p1_data_container.hpp"

#include <dune/grid/common/geometry.hh>

namespace Bempp {

    struct P1GridFamily {};

    class P1Grid {

        public: 

            typedef double ctype;

            enum {dimension = 2};
            enum {dimensionworld = 3};
    
    };

    

}


#endif
