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
#include "triangle_imp/bempp_grid_triangle_factory.hpp"

namespace BemppGrid {

    inline TriangleGrid::TriangleGrid(const shared_ptr<DataContainer>& data) : m_data(data) {

        for (int i = 0; i < data->levels(); ++i) {

            m_levelIndexSet.push_back(
                    shared_ptr<typename GridFamily::Traits::LevelIndexSet>(
                        new LevelIndexSetImp(data, i)));

        }


    }

    template <int cd, Dune::PartitionIteratorType pitype>
    typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator TriangleGrid::lbegin(int level) const {

        return typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator(LevelIteratorImp<cd, pitype, const TriangleGrid>(m_data, level, 0));

    }

    template <int cd, Dune::PartitionIteratorType pitype>
    typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator TriangleGrid::lend(int level) const {

        return typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator(LevelIteratorImp<cd, pitype, const TriangleGrid>(m_data, level, 
                    m_data->numberOfEntities<cd>(level)));

    }

    template <int cd, Dune::PartitionIteratorType pitype>
    typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator TriangleGrid::leafbegin() const {

        return typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator(LevelIteratorImp<cd, pitype, const TriangleGrid>(
                    m_data, m_data->levels() - 1, 0));

    }

    template <int cd, Dune::PartitionIteratorType pitype>
    typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator TriangleGrid::leafend() const {

        int level = m_data->levels() - 1;
        return typename TriangleGrid::Codim<cd>::template Partition<pitype>::LevelIterator(LevelIteratorImp<cd, pitype, const TriangleGrid>(m_data, level, 
                    m_data->numberOfEntities<cd>(level)));

    }

    inline const typename TriangleGrid::GridFamily::Traits::LevelIndexSet& TriangleGrid::levelIndexSet(int level) const {

        return *m_levelIndexSet[level];

    }

    inline const typename TriangleGrid::GridFamily::Traits::LeafIndexSet& TriangleGrid::leafIndexSet() const {

        return *m_levelIndexSet[m_data->levels() - 1];

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

    inline int TriangleGrid::maxLevel() const {

        return m_data->levels() - 1;

    }

    inline int TriangleGrid::size(int level, int codim) const {

        if (codim == 0)
            return m_data->numberOfEntities<0>(level);
        if (codim == 1)
            return m_data->numberOfEntities<1>(level);
        if (codim ==2)
            return m_data->numberOfEntities<2>(level);
        throw std::runtime_error("TriangleGrid::size(): Need 0 <= codim <= 2");


    }

    inline int TriangleGrid::size(int codim) const {

        return size(maxLevel(), codim);


    }

    inline int TriangleGrid::size(int level, Dune::GeometryType type) const {

        if (type.isVertex())
            return size(level, 2);
        if (type.isLine())
            return size(level, 1);
        if (type.isTriangle())
            return size(level, 0);
        throw std::runtime_error("TriangleGrid::size(): unknwon type");


    }

    inline int TriangleGrid::size(Dune::GeometryType type) const {

        return size(maxLevel(), type);

    }

}




#endif
