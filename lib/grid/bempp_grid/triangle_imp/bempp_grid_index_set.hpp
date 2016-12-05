#ifndef bempp_grid_triangle_imp_index_set_hpp
#define bempp_grid_triangle_imp_index_set_hpp

#include "../../../common/common.hpp"
#include <dune/grid/common/entity.hh>

namespace BemppGrid {

    class TriangleGrid;
    template<int, int, class> class EntityImp;

    class LevelIndexSetImp : public Dune::IndexSetDefaultImplementation<const TriangleGrid, LevelIndexSetImp>
    {

        public:

        typedef unsigned int IndexType;

        template <int cd>
        IndexType index(const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity& entity) const
        {

            return TriangleGrid::entityIndex<cd>(entity);

        }

        IndexType subIndex(const typename TriangleGrid::GridFamily::Traits::Codim<0>::Entity& entity,
                int i, unsigned int codim) const{

            if (codim == 0)
                return index<0>(entity);
            if (codim == 1)
                return index<1>(*entity.template subEntity<1>(i));
            if (codim == 2)
                return index<2>(*entity.template subEntity<2>(i));
            throw std::runtime_error("LevelIndexSetImp::subIndex(): Error. Must have 0 <= codim <= 2.");


        }



    }; 

}


#endif
