#ifndef bempp_triangle_grid_imp_hpp
#define bempp_triangle_grid_imp_hpp

#include "bempp_triangle_grid.hpp"
#include "triangle_imp/bempp_grid_data_container.hpp"
#include "triangle_imp/bempp_grid_entity_pointer.hpp"
#include "triangle_imp/bempp_grid_entity_seed.hpp"
#include "triangle_imp/bempp_grid_geometry.hpp"
#include "triangle_imp/bempp_grid_entity.hpp"
#include "triangle_imp/bempp_grid_entity_pointer.hpp"
#include "triangle_imp/bempp_grid_level_iterator.hpp"
#include "triangle_imp/bempp_grid_index_set.hpp"
#include "triangle_imp/bempp_grid_id_set.hpp"
#include "bempp_triangle_grid_factory.hpp"

namespace BemppGrid {

    TriangleGrid::TriangleGrid(const shared_ptr<DataContainer>& data) : m_data(data) {}

    template <int cd>
    typename TriangleGrid::Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator TriangleGrid::lbegin(int level) const {

        return typename TriangleGrid::Codim<cd>::Iterator(LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>(m_data, level, 0));

    }

    template <int cd>
    typename TriangleGrid::Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator TriangleGrid::lend(int level) const {

        return typename TriangleGrid::Codim<cd>::Iterator(LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>(m_data, level, 
                    m_data->numberOfEntities<cd>(level)));

    }

    template <int cd>
    typename TriangleGrid::Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator TriangleGrid::leafbegin() const {

        return typename TriangleGrid::Codim<cd>::Iterator(LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>(
                    m_data, m_data->levels() - 1, 0));

    }

    template <int cd>
    typename TriangleGrid::Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator TriangleGrid::leafend() const {

        int level = m_data->levels() - 1;
        return typename TriangleGrid::Codim<cd>::Iterator(LevelIteratorImp<cd, Dune::All_Partition, const TriangleGrid>(m_data, level, 
                    m_data->numberOfEntities<cd>(level)));

    }

    inline const typename TriangleGrid::GridFamily::Traits::LevelIndexSet& TriangleGrid::levelIndexSet(int level) const {

        return m_levelIndexSet;

    }

    inline const typename TriangleGrid::GridFamily::Traits::LeafIndexSet& TriangleGrid::leafIndexSet() const {

        return m_leafIndexSet;

    }

    inline const typename TriangleGrid::GridFamily::Traits::GlobalIdSet& TriangleGrid::globalIdSet() const {

        return m_globalIdSet;

    }

    inline const typename TriangleGrid::GridFamily::Traits::LocalIdSet& TriangleGrid::localIdSet() const {

        return m_localIdSet;

    }
    
    template <int cd>
    TriangleGrid::GridFamily::Traits::LevelIndexSet::IndexType TriangleGrid::entityIndex(
            const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity& entity) {

        return TriangleGrid::getRealImplementation(entity).m_index;

    }

    template <int cd>
    int TriangleGrid::entityLevel(
            const typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity& entity) {

        return TriangleGrid::getRealImplementation(entity).m_level;

    }

    template <int cd>
    const EntityImp<cd, 2, const TriangleGrid>& TriangleGrid::getEntityImp(
            const typename GridFamily::Traits::Codim<cd>::Entity& entity){

        return TriangleGrid::getRealImplementation(entity);

    }

}




#endif
