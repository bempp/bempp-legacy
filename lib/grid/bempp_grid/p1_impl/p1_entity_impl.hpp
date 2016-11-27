#ifndef bempp_p1_entity_impl_hpp
#define bempp_p1_entity_impl_hpp

#include "p1_entity.hpp"
#include "p1_entity_pointer.hpp"
#include "p1_data_container.hpp"
#include <dune/geometry/type.hh>


namespace BemppGrid {

    P1EntityImp<0, 2, const P1Grid>::P1EntityImp(const shared_ptr<P1DataContainer>& data,
        int level, std::size_t index) : m_data(data), m_level(level), m_index(index)
    {}

    P1EntityImp<0, 2, const P1Grid>::Geometry P1EntityImp<0, 2, const P1Grid>::geometry() const {

               return m_data->geometry(Entity(*this)); 
    }

    Dune::PartitionType P1EntityImp<0, 2, const P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    int P1EntityImp<0, 2, const P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1EntityImp<0, 2, const P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        return geometryType;

    }

    bool P1EntityImp<0, 2, const P1Grid>::equals(const P1EntityImp<0, 2, const P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    P1EntityImp<0, 2, const P1Grid>::EntitySeed P1EntityImp<0, 2, const P1Grid>::seed() const {

        return P1EntityImp<0, 2, const P1Grid>::EntitySeed(P1EntitySeedImp<0, const P1Grid>(m_level, m_index));

    }

    template<>
    P1EntityImp<0, 2, const P1Grid>::EntityPointer<1> P1EntityImp<0, 2, const P1Grid>::subEntity<1>(int i) const {

        return EntityPointer<1>(
                P1EntityImp<1, 2, const P1Grid>(m_data, m_level,
                    m_data->element2Edges(m_level, m_index)[i]));

    }

    template<>
    P1EntityImp<0, 2, const P1Grid>::EntityPointer<2> P1EntityImp<0, 2, const P1Grid>::subEntity<2>(int i) const {

        return EntityPointer<2>(
                P1EntityImp<2, 2, const P1Grid>(m_data, m_level,
                    m_data->elements(m_level)[m_index][i]));

    }

    template<>
    P1EntityImp<0, 2, const P1Grid>::EntityPointer<0> P1EntityImp<0, 2, const P1Grid>::subEntity<0>(int i) const {

        return EntityPointer<0>(*this);

    }


    P1EntityImp<0, 2, const P1Grid>::EntityPointer<0> P1EntityImp<0, 2, const P1Grid>::father() const {

        int fatherIndex = m_data->getElementFatherIndex(m_level, m_index);
        return EntityPointer<0>(
                P1EntityImp<0, 2, const P1Grid>(m_data, m_level, fatherIndex));

    }
    
    template<>
    int P1EntityImp<0, 2, const P1Grid>::count<0>() const {

        return 1;

    }

    template<>
    int P1EntityImp<0, 2, const P1Grid>::count<1>() const {

        return 3;

    }

    template<>
    int P1EntityImp<0, 2, const P1Grid>::count<2>() const {

        return 3;

    }

    bool P1EntityImp<0, 2, const P1Grid>::hasFather() const {

        return (m_level > 0);

    } 

    bool P1EntityImp<0, 2, const P1Grid>::isLeaf() const {

        return (m_level == m_data->levels() - 1);

    } 

    bool P1EntityImp<0, 2, const P1Grid>::isRegular() const {

        return true;

    } 

    bool P1EntityImp<0, 2, const P1Grid>::isNew() const {

        if (m_level == 0) return false;

        int fatherIndex = m_data->getElementFatherIndex(m_level, m_index);
        return m_data->getElementSons(m_level - 1, fatherIndex).size() > 0;

    }

    P1EntityImp<1, 2, const P1Grid>::P1EntityImp(const shared_ptr<P1DataContainer>& data,
        int level, std::size_t index) : m_data(data), m_level(level), m_index(index)
    {}


    P1EntityImp<1, 2, const P1Grid>::Geometry P1EntityImp<1, 2, const P1Grid>::geometry() const {

               return m_data->geometry(Entity(*this)); 
    }

    int P1EntityImp<1, 2, const P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1EntityImp<1, 2, const P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeLine();
        return geometryType;

    }

    Dune::PartitionType P1EntityImp<1, 2, const P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool P1EntityImp<1, 2, const P1Grid>::equals(const P1EntityImp<1, 2, const P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    P1EntityImp<1, 2, const P1Grid>::EntitySeed P1EntityImp<1, 2, const P1Grid>::seed() const {

        return P1EntityImp<1, 2, const P1Grid>::EntitySeed(P1EntitySeedImp<1, const P1Grid>(m_level, m_index));

    }

    P1EntityImp<2, 2, const P1Grid>::P1EntityImp(const shared_ptr<P1DataContainer>& data,
        int level, std::size_t index) : m_data(data), m_level(level), m_index(index)
    {
    }

    P1EntityImp<2, 2, const P1Grid>::EntitySeed P1EntityImp<2, 2, const P1Grid>::seed() const {

        return P1EntityImp<2, 2, const P1Grid>::EntitySeed(P1EntitySeedImp<2, const P1Grid>(m_level, m_index));

    }

    P1EntityImp<2, 2, const P1Grid>::Geometry P1EntityImp<2, 2, const P1Grid>::geometry() const {

               return m_data->geometry(Entity(*this)); 
    }

    int P1EntityImp<2, 2, const P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1EntityImp<2, 2, const P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeVertex();
        return geometryType;

    }

    Dune::PartitionType P1EntityImp<2, 2, const P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool P1EntityImp<2, 2, const P1Grid>::equals(const P1EntityImp<2, 2, const P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }
}



#endif
