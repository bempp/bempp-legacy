#ifndef bempp_p1_level_iterator_hpp
#define bempp_p1_level_iterator_hpp

#include "p1_entity.hpp"
#include "p1_entity_pointer.hpp"
#include <dune/grid/common/gridenums.hh>

namespace BemppGrid {

class P1DataContainer;    
class P1Grid;

template<int, int, class> class P1EntityImp;

template <int codim, Dune::PartitionIteratorType, class> 
class P1LevelIteratorImp : public P1EntityPointerImp<codim, P1Grid> {

    public:

        typedef Dune::Entity<codim, 2, P1Grid, P1EntityImp> Entity;
        typedef P1EntityPointerImp<codim, P1Grid> Base;

        P1LevelIteratorImp(const shared_ptr<P1DataContainer>& data, int level, std::size_t index) :
            m_data(data), m_level(level), m_index(index), 
            Base(P1EntityPointerImp<codim, P1Grid>(P1EntityImp<codim, 2, P1Grid>(data, level, index))) {};


        void increment() const {

            m_index++;
            static_cast<const Base*>(this)->setEntity(
                    P1EntityImp<codim, 2, P1Grid>(m_data, m_level, m_index));

        }

            

    private:

        shared_ptr<P1DataContainer> m_data;
        int m_level;
        mutable std::size_t m_index;



};

}




#endif
