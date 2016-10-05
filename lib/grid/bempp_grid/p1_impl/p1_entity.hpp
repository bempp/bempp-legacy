#ifndef bempp_p1_entity_hpp
#define bempp_p1_entity_hpp

#include "../../../common/common.hpp"
#include "../bempp_grid_types.hpp"
#include "p1_geometry.hpp"
#include "p1_entity_seed.hpp"
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/common/entityseed.hh>
#include <dune/grid/common/entitypointer.hh>

namespace BemppGrid {

    class P1Grid;
    class P1DataContainer;

    template <int cd, int dim, class GridImp>
    class P1EntityImp {};

    template <int codim>
    class P1EntityPointerImp;

    // Elements
    template<>
    class P1EntityImp<0, 2, P1Grid> : 
    public Dune::EntityDefaultImplementation<0, 2, P1Grid, P1EntityImp>  {

        friend class P1EntitySeedImp<0>;

        public:

            typedef Dune::EntitySeed<P1Grid, P1EntitySeedImp<0>> EntitySeed;
            template <int codim> using 
                EntityPointer = Dune::EntityPointer<P1Grid, P1EntityPointerImp<codim>>;

            P1EntityImp(const shared_ptr<P1DataContainer>& data, 
                    int level, std::size_t index);

            int level() const;
            Dune::Geometry<2, 3, P1Grid, P1GridGeometry> geometry() const;
            Dune::GeometryType type() const;
            Dune::PartitionType partitionType() const;
            bool equals(const P1EntityImp<0, 2, P1Grid>& other) const;
            EntitySeed seed() const; 

            template<int cd>
            EntityPointer<cd> subEntity(int i) const;


        private:

            shared_ptr<P1DataContainer> m_data;
            int m_level;
            std::size_t m_index;
            P1GridGeometry<2, 3, P1Grid> m_geometry;

    };

    // Edges
    template<>
    class P1EntityImp<1, 2, P1Grid> : 
    public Dune::EntityDefaultImplementation<1, 2, P1Grid, P1EntityImp> {

        public:
            typedef Dune::EntitySeed<P1Grid, P1EntitySeedImp<1>> EntitySeed;

            P1EntityImp(const shared_ptr<P1DataContainer>& data,
                    int level, std::size_t index);

            int level() const;
            Dune::Geometry<1, 3, P1Grid, P1GridGeometry> geometry() const;
            Dune::GeometryType type() const;
            Dune::PartitionType partitionType() const;
            bool equals(const P1EntityImp<1, 2, P1Grid>& other) const;
            EntitySeed seed() const;

        private:

            shared_ptr<P1DataContainer> m_data;
            int m_level;
            std::size_t m_index;
            P1GridGeometry<1, 3, P1Grid> m_geometry;

    };

    // Vertices
    template<>
    class P1EntityImp<2, 2, P1Grid> :
    public Dune::EntityDefaultImplementation<2, 2, P1Grid, P1EntityImp> {

        public:
            typedef Dune::EntitySeed<P1Grid, P1EntitySeedImp<2>> EntitySeed;

            P1EntityImp(const shared_ptr<P1DataContainer>& data,
                    int level, std::size_t index);

            int level() const;
            Dune::Geometry<0, 3, P1Grid, P1GridGeometry> geometry() const;
            Dune::GeometryType type() const;
            Dune::PartitionType partitionType() const;
            bool equals(const P1EntityImp<2, 2, P1Grid>& other) const;
            EntitySeed seed() const;

        private:

            shared_ptr<P1DataContainer> m_data;
            int m_level;
            std::size_t m_index;
            P1GridGeometry<0, 3, P1Grid> m_geometry;
    };

}

#include "p1_entity_impl.hpp"

#endif
