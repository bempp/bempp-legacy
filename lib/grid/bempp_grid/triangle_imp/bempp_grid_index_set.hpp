#ifndef bempp_grid_triangle_imp_index_set_hpp
#define bempp_grid_triangle_imp_index_set_hpp

#include "../../../common/common.hpp"
#include "bempp_grid_data_container.hpp"
#include <dune/grid/common/entity.hh>
#include <dune/geometry/type.hh>

namespace BemppGrid {

    class TriangleGrid;
    class DataContainer;
    template<int, int, class> class EntityImp;

    class LevelIndexSetImp : public Dune::IndexSetDefaultImplementation<const TriangleGrid, LevelIndexSetImp>
    {

        public:

        typedef unsigned int IndexType;

        LevelIndexSetImp(const shared_ptr<DataContainer>& data, unsigned int level) : m_data(data), m_level(level) {}

        template <int cd>
        IndexType index(const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity& entity) const
        {

            return TriangleGrid::entityIndex<cd>(entity);

        }

        template <int cd>
        IndexType subIndex(const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity& entity,
                int i, unsigned int codim) const{

            if (codim > cd) 
                throw std::runtime_error("LevelIndexSetImp::subIndex(): Error Must have 0 <=codim <= cd.");

            return index<cd>(*entity.template subEntity<cd>(i));

        }

        IndexType size(Dune::GeometryType type) const {

            if (type.isVertex())
                return m_data->numberOfEntities<2>(m_level);
            if (type.isLine())
                return m_data->numberOfEntities<1>(m_level);
            if (type.isTriangle())
                return m_data->numberOfEntities<0>(m_level);

            throw std::runtime_error("LevelIndexSetImp::size(): Unknownn Geometry type");


        }

        IndexType size(int codim) const {

            if (codim == 2)
                return m_data->numberOfEntities<2>(m_level);
            if (codim == 1)
                return m_data->numberOfEntities<1>(m_level);
            if (codim == 0)
                return m_data->numberOfEntities<0>(m_level);

            throw std::runtime_error("LevelIndexSetImp::size(): Unknownn Geometry type");


        }

        template <class EntityType>
        bool contains(const EntityType& e) const {

            return e.level() == m_level;

        }

        private:

        shared_ptr<DataContainer> m_data;
        unsigned int m_level;

    }; 

}


#endif
