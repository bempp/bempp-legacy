#ifndef bempp_p1_index_set_hpp
#define bempp_p1_index_set_hpp

#include "../../../common/common.hpp"
#include <dune/grid/common/entity.hh>

namespace BemppGrid {

    class P1Grid;
    template<int, int, class> class P1EntityImp;

    class P1LevelIndexSetImp : public Dune::IndexSetDefaultImplementation<P1Grid, P1LevelIndexSetImp>
    {

        public:

        typedef std::size_t IndexType;

        template <int cd>
        IndexType index(const Dune::Entity<cd, 2, const P1Grid,P1EntityImp>& entity)
        {

            return P1Grid::entityIndex<cd>(entity);

        }





    }; 

}


#endif
