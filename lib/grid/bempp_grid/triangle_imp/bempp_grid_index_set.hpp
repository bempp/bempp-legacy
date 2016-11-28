#ifndef bempp_grid_triangle_imp_index_set_hpp
#define bempp_grid_triangle_imp_index_set_hpp

#include "../../../common/common.hpp"
#include <dune/grid/common/entity.hh>

namespace BemppGrid {

    class TriangleGrid;
    template<int, int, class> class EntityImp;

    class LevelIndexSetImp : public Dune::IndexSetDefaultImplementation<TriangleGrid, LevelIndexSetImp>
    {

        public:

        typedef std::size_t IndexType;

        template <int cd>
        IndexType index(const Dune::Entity<cd, 2, const TriangleGrid,EntityImp>& entity)
        {

            return TriangleGrid::entityIndex<cd>(entity);

        }





    }; 

}


#endif
