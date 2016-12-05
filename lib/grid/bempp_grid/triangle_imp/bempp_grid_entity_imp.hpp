#ifndef bempp_grid_triangle_imp_entity_imp_hpp
#define bempp_grid_triangle_imp_entity_imp_hpp

#include "bempp_grid_entity.hpp"
#include "bempp_grid_data_container.hpp"
#include <dune/geometry/type.hh>


namespace BemppGrid {

    EntityImp<0, 2, const TriangleGrid>::EntityImp(const shared_ptr<DataContainer>& data,
        int level, unsigned int index) : m_data(data), m_level(level), m_index(index)
    {}

    EntityImp<0, 2, const TriangleGrid>::Geometry EntityImp<0, 2, const TriangleGrid>::geometry() const {

               return m_data->geometry<0>(Entity(*this)); 
    }

    Dune::PartitionType EntityImp<0, 2, const TriangleGrid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    int EntityImp<0, 2, const TriangleGrid>::level() const {
        return m_level;
    }     

    Dune::GeometryType EntityImp<0, 2, const TriangleGrid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        return geometryType;

    }

    bool EntityImp<0, 2, const TriangleGrid>::equals(const EntityImp<0, 2, const TriangleGrid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    EntityImp<0, 2, const TriangleGrid>::EntitySeed EntityImp<0, 2, const TriangleGrid>::seed() const {

        return EntityImp<0, 2, const TriangleGrid>::EntitySeed(EntitySeedImp<0, const TriangleGrid>(m_level, m_index));

    }

    template<>
    EntityImp<0, 2, const TriangleGrid>::EntityPointer<1> EntityImp<0, 2, const TriangleGrid>::subEntity<1>(int i) const {

        return EntityPointer<1>(
                EntityImp<1, 2, const TriangleGrid>(m_data, m_level,
                    m_data->element2Edges(m_level, m_index)[i]));

    }

    template<>
    EntityImp<0, 2, const TriangleGrid>::EntityPointer<2> EntityImp<0, 2, const TriangleGrid>::subEntity<2>(int i) const {

        return EntityPointer<2>(
                EntityImp<2, 2, const TriangleGrid>(m_data, m_level,
                    m_data->elements(m_level)[m_index][i]));

    }

    template<>
    EntityImp<0, 2, const TriangleGrid>::EntityPointer<0> EntityImp<0, 2, const TriangleGrid>::subEntity<0>(int i) const {

        return EntityPointer<0>(*this);

    }


    EntityImp<0, 2, const TriangleGrid>::EntityPointer<0> EntityImp<0, 2, const TriangleGrid>::father() const {

        int fatherIndex = m_data->getElementFatherIndex(m_level, m_index);
        return EntityPointer<0>(
                EntityImp<0, 2, const TriangleGrid>(m_data, m_level, fatherIndex));

    }
    
    template<>
    int EntityImp<0, 2, const TriangleGrid>::count<0>() const {

        return 1;

    }

    template<>
    int EntityImp<0, 2, const TriangleGrid>::count<1>() const {

        return 3;

    }

    template<>
    int EntityImp<0, 2, const TriangleGrid>::count<2>() const {

        return 3;

    }

    bool EntityImp<0, 2, const TriangleGrid>::hasFather() const {

        return (m_level > 0);

    } 

    bool EntityImp<0, 2, const TriangleGrid>::isLeaf() const {

        return (m_level == m_data->levels() - 1);

    } 

    bool EntityImp<0, 2, const TriangleGrid>::isRegular() const {

        return true;

    } 

    bool EntityImp<0, 2, const TriangleGrid>::isNew() const {

        if (m_level == 0) return false;

        int fatherIndex = m_data->getElementFatherIndex(m_level, m_index);
        return m_data->getElementSons(m_level - 1, fatherIndex).size() > 0;

    }


    EntityImp<1, 2, const TriangleGrid>::EntityImp(const shared_ptr<DataContainer>& data,
        int level, unsigned int index) : m_data(data), m_level(level), m_index(index)
    {}


    EntityImp<1, 2, const TriangleGrid>::Geometry EntityImp<1, 2, const TriangleGrid>::geometry() const {

               return m_data->geometry<1>(Entity(*this)); 
    }

    int EntityImp<1, 2, const TriangleGrid>::level() const {
        return m_level;
    }     

    Dune::GeometryType EntityImp<1, 2, const TriangleGrid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeLine();
        return geometryType;

    }

    Dune::PartitionType EntityImp<1, 2, const TriangleGrid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool EntityImp<1, 2, const TriangleGrid>::equals(const EntityImp<1, 2, const TriangleGrid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    EntityImp<1, 2, const TriangleGrid>::EntitySeed EntityImp<1, 2, const TriangleGrid>::seed() const {

        return EntityImp<1, 2, const TriangleGrid>::EntitySeed(EntitySeedImp<1, const TriangleGrid>(m_level, m_index));

    }

    EntityImp<2, 2, const TriangleGrid>::EntityImp(const shared_ptr<DataContainer>& data,
        int level, unsigned int index) : m_data(data), m_level(level), m_index(index)
    {
    }

    EntityImp<2, 2, const TriangleGrid>::EntitySeed EntityImp<2, 2, const TriangleGrid>::seed() const {

        return EntityImp<2, 2, const TriangleGrid>::EntitySeed(EntitySeedImp<2, const TriangleGrid>(m_level, m_index));

    }

    EntityImp<2, 2, const TriangleGrid>::Geometry EntityImp<2, 2, const TriangleGrid>::geometry() const {

               return m_data->geometry<2>(Entity(*this)); 
    }

    int EntityImp<2, 2, const TriangleGrid>::level() const {
        return m_level;
    }     

    Dune::GeometryType EntityImp<2, 2, const TriangleGrid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeVertex();
        return geometryType;

    }

    Dune::PartitionType EntityImp<2, 2, const TriangleGrid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool EntityImp<2, 2, const TriangleGrid>::equals(const EntityImp<2, 2, const TriangleGrid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    unsigned int EntityImp<0, 2, const TriangleGrid>::id() const {

        return m_data->id<0>(Entity(*this));

    }

    unsigned int EntityImp<1, 2, const TriangleGrid>::id() const {

        return m_data->id<1>(Entity(*this));

    }

    unsigned int EntityImp<2, 2, const TriangleGrid>::id() const {

        return m_data->id<2>(Entity(*this));

    }
}



#endif
