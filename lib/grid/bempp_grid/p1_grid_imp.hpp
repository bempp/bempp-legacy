#ifndef bempp_p1_grid_imp_hpp
#define bempp_p1_grid_imp_hpp

#include "p1_grid.hpp"
#include "p1_impl/p1_data_container.hpp"
#include "p1_impl/p1_entity_seed.hpp"
#include "p1_impl/p1_geometry.hpp"
#include "p1_impl/p1_entity.hpp"
#include "p1_impl/p1_entity_pointer.hpp"
#include "p1_impl/p1_data_container.hpp"
#include "p1_impl/p1_level_iterator.hpp"
#include "p1_impl/p1_index_set.hpp"

namespace BemppGrid {

    P1Grid::P1Grid(const shared_ptr<P1DataContainer>& data) : m_data(data) {}

    template <int cd>
    typename P1Grid::Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator P1Grid::lbegin(int level) const {

        return typename P1Grid::Codim<cd>::Iterator(P1LevelIteratorImp<cd, Dune::All_Partition, const P1Grid>(m_data, level, 0));

    }

    template <int cd>
    typename P1Grid::Codim<cd>::template Partition<Dune::All_Partition>::LevelIterator P1Grid::lend(int level) const {

        return typename P1Grid::Codim<cd>::Iterator(P1LevelIteratorImp<cd, Dune::All_Partition, const P1Grid>(m_data, level, 
                    m_data->numberOfEntities<cd>(level)));

    }
    
    template <int cd>
    P1Grid::GridFamily::Traits::LevelIndexSet::IndexType P1Grid::entityIndex(
            const typename P1Grid::GridFamily::Traits::Codim<cd>::Entity& entity) {

        return P1Grid::getRealImplementation(entity).m_index;

    }

    template <int cd>
    int P1Grid::entityLevel(
            const typename P1Grid::GridFamily::Traits::Codim<cd>::Entity& entity) {

        return P1Grid::getRealImplementation(entity).m_level;

    }

}




#endif
