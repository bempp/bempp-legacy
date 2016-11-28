#ifndef bempp_grid_triangle_imp_id_set_hpp
#define bempp_grid_triangle_imp_id_set_hpp

#include "../../../common/common.hpp"
#include <dune/grid/common/entity.hh>

namespace BemppGrid {

    class TriangleGrid;
    template<int, int, class> class EntityImp;

    class IdSetImp : public Dune::IdSet<const TriangleGrid, IdSetImp, std::size_t> 
    {

        public:

        typedef std::size_t IdType;

        template <int cd>
        IdType id(const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity& entity) const
        {

            return TriangleGrid::getEntityImp<cd>(entity).id();

        }

        inline IdType subId(const typename TriangleGrid::GridFamily::Traits::Codim<0>::Entity& entity,
                int i, unsigned int codim) const {

            if (codim == 0)
                return id<0>(entity);
            if (codim == 1)
                return id<1>(*entity.template subEntity<1>(i));
            if (codim == 2)
                return id<2>(*entity.template subEntity<2>(i));
            throw std::runtime_error("IdSetImp::subId(): Error. Must have 0 <= codim <= 2.");


        }


    }; 

}


#endif
