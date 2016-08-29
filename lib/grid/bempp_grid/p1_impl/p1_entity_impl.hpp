#ifndef bempp_p1_entity_impl_hpp
#define bempp_p1_entity_impl_hpp

#include "p1_entity.hpp"
#include "p1_data_container.hpp"
#include <dune/geometry/type.hh>


namespace BemppGrid {

    P1Entity<0, 2, P1Grid>::P1Entity(const shared_ptr<P1DataContainer>& data,
        int level, int index) : m_data(data), m_level(level), m_index(index)
    {}

    Dune::Geometry<2, 3, P1Grid, P1GridGeometry> P1Entity<0, 2, P1Grid>::geometry() const {
        const auto& nodeIndices = m_data->elements(m_level)[m_index]; 
        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        std::vector<Dune::FieldVector<double, 3>> vertices;
        for (int i = 0; i < 3; ++i) vertices.push_back(m_data->nodes(m_level)[nodeIndices[i]]);
        return Dune::Geometry<2, 3, P1Grid, P1GridGeometry>(
                P1GridGeometry<2, 3, P1Grid>(geometryType, vertices));

    }

    Dune::PartitionType P1Entity<0, 2, P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    int P1Entity<0, 2, P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1Entity<0, 2, P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        return geometryType;

    }

    bool P1Entity<0, 2, P1Grid>::equals(const P1Entity<0, 2, P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    P1Entity<1, 2, P1Grid>::P1Entity(const shared_ptr<P1DataContainer>& data,
        int level, int index) : m_data(data), m_level(level), m_index(index)
    {}


    Dune::Geometry<1, 3, P1Grid, P1GridGeometry> P1Entity<1, 2, P1Grid>::geometry() const {
        const auto& nodeIndices = m_data->elements(m_level)[m_index]; 
        Dune::GeometryType geometryType;
        geometryType.makeLine();
        std::vector<Dune::FieldVector<double, 3>> vertices;
        for (int i = 0; i < 3; ++i) vertices.push_back(m_data->nodes(m_level)[nodeIndices[i]]);
        return Dune::Geometry<1, 3, P1Grid, P1GridGeometry>(
                P1GridGeometry<1, 3, P1Grid>(geometryType, vertices));

    }

    int P1Entity<1, 2, P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1Entity<1, 2, P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeLine();
        return geometryType;

    }

    Dune::PartitionType P1Entity<1, 2, P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool P1Entity<1, 2, P1Grid>::equals(const P1Entity<1, 2, P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    P1Entity<2, 2, P1Grid>::P1Entity(const shared_ptr<P1DataContainer>& data,
        int level, int index) : m_data(data), m_level(level), m_index(index)
    {}


    Dune::Geometry<0, 3, P1Grid, P1GridGeometry> P1Entity<2, 2, P1Grid>::geometry() const {
        const auto& nodeIndices = m_data->elements(m_level)[m_index]; 
        Dune::GeometryType geometryType;
        geometryType.makeVertex();
        std::vector<Dune::FieldVector<double, 3>> vertices;
        for (int i = 0; i < 3; ++i) vertices.push_back(m_data->nodes(m_level)[nodeIndices[i]]);
        return Dune::Geometry<0, 3, P1Grid, P1GridGeometry>(
                P1GridGeometry<0, 3, P1Grid>(geometryType, vertices));

    }

    int P1Entity<2, 2, P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1Entity<2, 2, P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeVertex();
        return geometryType;

    }

    Dune::PartitionType P1Entity<2, 2, P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool P1Entity<2, 2, P1Grid>::equals(const P1Entity<2, 2, P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }
}



#endif
