#ifndef bempp_p1_entity_impl_hpp
#define bempp_p1_entity_impl_hpp

#include "p1_entity.hpp"
#include "p1_entity_pointer.hpp"
#include "p1_data_container.hpp"
#include <dune/geometry/type.hh>


namespace BemppGrid {

    P1EntityImp<0, 2, P1Grid>::P1EntityImp(const shared_ptr<P1DataContainer>& data,
        int level, std::size_t index) : m_data(data), m_level(level), m_index(index)
    {}

    Dune::Geometry<2, 3, P1Grid, P1GridGeometry> P1EntityImp<0, 2, P1Grid>::geometry() const {

        if (m_geometry == nullptr) {
            const auto& nodeIndices = m_data->elements(m_level)[m_index]; 
            Dune::GeometryType geometryType;
            geometryType.makeTriangle();
            std::vector<Dune::FieldVector<double, 3>> vertices;
            for (int i = 0; i < 3; ++i) vertices.push_back(m_data->nodes(m_level)[nodeIndices[i]]);
            m_geometry = shared_ptr<P1GridGeometry<2, 3, P1Grid>>(
                    new P1GridGeometry<2, 3, P1Grid>(geometryType, vertices));
        }

        return Dune::Geometry<2, 3, P1Grid, P1GridGeometry>(*m_geometry);
                
    }

    Dune::PartitionType P1EntityImp<0, 2, P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    int P1EntityImp<0, 2, P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1EntityImp<0, 2, P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        return geometryType;

    }

    bool P1EntityImp<0, 2, P1Grid>::equals(const P1EntityImp<0, 2, P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    P1EntityImp<0, 2, P1Grid>::EntitySeed P1EntityImp<0, 2, P1Grid>::seed() const {

        return P1EntityImp<0, 2, P1Grid>::EntitySeed(P1EntitySeedImp<0>(m_level, m_index));

    }

    template<>
    P1EntityImp<0, 2, P1Grid>::EntityPointer<1> P1EntityImp<0, 2, P1Grid>::subEntity<1>(int i) const {

        return EntityPointer<1>(
                P1EntityImp<1, 2, P1Grid>(m_data, m_level,
                    m_data->element2Edges(m_level, m_index)[i]));

    }

    template<>
    P1EntityImp<0, 2, P1Grid>::EntityPointer<2> P1EntityImp<0, 2, P1Grid>::subEntity<2>(int i) const {

        return EntityPointer<2>(
                P1EntityImp<2, 2, P1Grid>(m_data, m_level,
                    m_data->elements(m_level)[m_index][i]));

    }

    template<>
    P1EntityImp<0, 2, P1Grid>::EntityPointer<0> P1EntityImp<0, 2, P1Grid>::subEntity<0>(int i) const {

        return EntityPointer<0>(*this);

    }
    

    P1EntityImp<1, 2, P1Grid>::P1EntityImp(const shared_ptr<P1DataContainer>& data,
        int level, std::size_t index) : m_data(data), m_level(level), m_index(index)
    {}


    Dune::Geometry<1, 3, P1Grid, P1GridGeometry> P1EntityImp<1, 2, P1Grid>::geometry() const {

        if (m_geometry == nullptr) {
            const auto& nodeIndices = m_data->edges(m_level)[m_index]; 
            Dune::GeometryType geometryType;
            geometryType.makeLine();
            std::vector<Dune::FieldVector<double, 3>> vertices;
            for (int i = 0; i < 2; ++i) vertices.push_back(m_data->nodes(m_level)[nodeIndices[i]]);
            m_geometry = shared_ptr<P1GridGeometry<1, 3, P1Grid>>(
                    new P1GridGeometry<1, 3, P1Grid>(geometryType, vertices));

        }
        return Dune::Geometry<1, 3, P1Grid, P1GridGeometry>(*m_geometry);

    }

    int P1EntityImp<1, 2, P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1EntityImp<1, 2, P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeLine();
        return geometryType;

    }

    Dune::PartitionType P1EntityImp<1, 2, P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool P1EntityImp<1, 2, P1Grid>::equals(const P1EntityImp<1, 2, P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }

    P1EntityImp<1, 2, P1Grid>::EntitySeed P1EntityImp<1, 2, P1Grid>::seed() const {

        return P1EntityImp<1, 2, P1Grid>::EntitySeed(P1EntitySeedImp<1>(m_level, m_index));

    }

    P1EntityImp<2, 2, P1Grid>::P1EntityImp(const shared_ptr<P1DataContainer>& data,
        int level, std::size_t index) : m_data(data), m_level(level), m_index(index)
    {
    }

    P1EntityImp<2, 2, P1Grid>::EntitySeed P1EntityImp<2, 2, P1Grid>::seed() const {

        return P1EntityImp<2, 2, P1Grid>::EntitySeed(P1EntitySeedImp<2>(m_level, m_index));

    }

    Dune::Geometry<0, 3, P1Grid, P1GridGeometry> P1EntityImp<2, 2, P1Grid>::geometry() const {

        if (m_geometry == nullptr) {
            const auto& node = m_data->nodes(m_level)[m_index]; 
            Dune::GeometryType geometryType;
            geometryType.makeVertex();
            std::vector<Dune::FieldVector<double, 3>> vertices;
            vertices.push_back(node);
            std::cout << "Output " << node << " " << vertices[0] << std::endl;
            m_geometry = shared_ptr<P1GridGeometry<0, 3, P1Grid>>(
                    new P1GridGeometry<0, 3, P1Grid>(geometryType, vertices));
        }
        return Dune::Geometry<0, 3, P1Grid, P1GridGeometry>(*m_geometry);

    }

    int P1EntityImp<2, 2, P1Grid>::level() const {
        return m_level;
    }     

    Dune::GeometryType P1EntityImp<2, 2, P1Grid>::type() const {

        Dune::GeometryType geometryType;
        geometryType.makeVertex();
        return geometryType;

    }

    Dune::PartitionType P1EntityImp<2, 2, P1Grid>::partitionType() const {

        return Dune::InteriorEntity;

    }

    bool P1EntityImp<2, 2, P1Grid>::equals(const P1EntityImp<2, 2, P1Grid>& other) const {

        return m_level == other.m_level && m_index == other.m_index;
    }
}



#endif
