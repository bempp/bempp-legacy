#ifndef bempp_p1_entity_hpp
#define bempp_p1_entity_hpp

#include "../../../common/common.hpp"
#include "../bempp_grid_types.hpp"
#include "p1_geometry.hpp"
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/entity.hh>


namespace BemppGrid {

    class P1Grid;
    class P1DataContainer;

    template <int cd, int dim, class GridImp>
    class P1EntityImp {};

    // Elements
    template<>
    class P1EntityImp<0, 2, P1Grid> : 
    public Dune::EntityDefaultImplementation<0, 2, P1Grid, P1EntityImp>  {

        public:

            P1EntityImp(const shared_ptr<P1DataContainer>& data, 
                    int level, int index);

            int level() const;
            Dune::Geometry<2, 3, P1Grid, P1GridGeometry> geometry() const;
            Dune::GeometryType type() const;
            Dune::PartitionType partitionType() const;
            bool equals(const P1EntityImp<0, 2, P1Grid>& other) const;

        private:

            shared_ptr<P1DataContainer> m_data;
            int m_level;
            int m_index;

    };

    // Edges
    template<>
    class P1EntityImp<1, 2, P1Grid> : 
    public Dune::EntityDefaultImplementation<1, 2, P1Grid, P1EntityImp> {

        public:

            P1EntityImp(const shared_ptr<P1DataContainer>& data,
                    int level, int index);

            int level() const;
            Dune::Geometry<1, 3, P1Grid, P1GridGeometry> geometry() const;
            Dune::GeometryType type() const;
            Dune::PartitionType partitionType() const;
            bool equals(const P1EntityImp<1, 2, P1Grid>& other) const;

        private:

            shared_ptr<P1DataContainer> m_data;
            int m_level;
            int m_index;

    };

    // Vertices
    template<>
    class P1EntityImp<2, 2, P1Grid> :
    public Dune::EntityDefaultImplementation<2, 2, P1Grid, P1EntityImp> {

        public:

            P1EntityImp(const shared_ptr<P1DataContainer>& data,
                    int level, int index);

            int level() const;
            Dune::Geometry<0, 3, P1Grid, P1GridGeometry> geometry() const;
            Dune::GeometryType type() const;
            Dune::PartitionType partitionType() const;
            bool equals(const P1EntityImp<2, 2, P1Grid>& other) const;

        private:

            shared_ptr<P1DataContainer> m_data;
            int m_level;
            int m_index;
    };

}

#include "p1_entity_impl.hpp"

#endif
