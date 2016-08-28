#ifndef bempp_p1_entity_impl_hpp
#define bempp_p1_entity_impl_hpp

#include "p1_entity.hpp"

namespace Bempp {

    P1Entity<0, 2, P1Grid>::P1Entity(const shared_ptr<P1DataContainer>& data,
        int level, int index) : m_data(data), m_level(level), m_index(index)
    {}

    int P1Entity<0, 2, P1Grid>::level() const {
        return m_level;
    }     

    P1Entity<1, 2, P1Grid>::P1Entity(const shared_ptr<P1DataContainer>& data,
        int level, int index) : m_data(data), m_level(level), m_index(index)
    {}

    int P1Entity<1, 2, P1Grid>::level() const {
        return m_level;
    }     

    P1Entity<2, 2, P1Grid>::P1Entity(const shared_ptr<P1DataContainer>& data,
        int level, int index) : m_data(data), m_level(level), m_index(index)
    {}

    int P1Entity<2, 2, P1Grid>::level() const {
        return m_level;
    }     
}



#endif
