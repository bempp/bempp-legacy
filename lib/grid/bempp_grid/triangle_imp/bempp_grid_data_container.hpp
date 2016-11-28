#ifndef bempp_grid_triangle_imp_data_container_hpp
#define bempp_grid_triangle_imp_data_container_hpp

#include "../../../common/common.hpp"
#include "../bempp_grid_types.hpp"
#include "bempp_grid_geometry.hpp"
#include <dune/common/fvector.hh>
#include <vector>
#include <array>
#include <tuple>

namespace BemppGrid {

    class TriangleGrid;

    class DataContainer {

        public:

            typedef std::vector<Dune::FieldVector<double, 3>> NodesContainer;
            typedef std::vector<std::array<std::size_t, 3>> ElementsContainer;
            typedef std::vector<std::array<std::size_t, 2>> EdgesContainer;

            template <int cd>
            using DuneGeometry = typename TriangleGrid::GridFamily::Traits::Codim<cd>::Geometry;

            template <int cd>
            using Geometry = BemppGrid::Geometry<2-cd, 3,const TriangleGrid>;

            template <int cd>
            using Entity = typename TriangleGrid::GridFamily::Traits::Codim<cd>::Entity;

            DataContainer();
            ~DataContainer();

            void init( const shared_ptr<NodesContainer>& nodes, 
                    const shared_ptr<ElementsContainer>& elements);

            const NodesContainer& nodes(int level) const;
            const ElementsContainer& elements(int level) const;
            const EdgesContainer& edges(int level) const;

            int numberOfNodes(int level) const;
            int numberOfElements(int level) const;
            int numberOfEdges(int level) const;
            int levels() const;

            template <int cd>
            int numberOfEntities(int level) const;

            template <int cd>
            DuneGeometry<cd> geometry(const Entity<cd>& entity);

            const std::array<std::size_t, 3>& element2Edges(int level, std::size_t elementIndex) const;
            const std::vector<size_t>& edge2Elements(int level, std::size_t edgeIndex) const;
            const std::vector<size_t>& node2Elements(int level, std::size_t nodeIndex) const;
            const std::vector<size_t>& node2Edges(int level, std::size_t nodeIndex) const; 

            const NodesContainer getEntityNodes(int codim, int level, int index) const;

            int getElementFatherIndex(int level, std::size_t elementIndex) const;
            const std::vector<size_t>& getElementSons(int level, std::size_t elementIndex) const;

        private:

            template <int cd>
            void computeGeometries(int level);

            int m_levels;
            std::vector<shared_ptr<NodesContainer>> m_nodes;
            std::vector<shared_ptr<ElementsContainer>> m_elements;
            std::vector<EdgesContainer> m_edges;

            std::vector<std::vector<std::size_t>> m_nodeIds;
            std::vector<std::vector<std::size_t>> m_elementIds;
            std::vector<std::vector<std::size_t>> m_edgeIds;

            std::size_t m_idCounter;

            std::vector<std::vector<std::array<std::size_t, 3>>> m_element2Edges;
            std::vector<std::vector<std::vector<std::size_t>>> m_edge2Elements;
            std::vector<std::vector<std::vector<std::size_t>>> m_node2Elements;
            std::vector<std::vector<std::vector<std::size_t>>> m_node2Edges;

            std::vector<std::vector<std::size_t>> m_fatherElements;
            std::vector<std::vector<std::vector<std::size_t>>> m_sonElements;

            std::tuple<
                std::vector<shared_ptr<std::vector<Geometry<0>>>>,
                std::vector<shared_ptr<std::vector<Geometry<1>>>>,
                std::vector<shared_ptr<std::vector<Geometry<2>>>>
                    > m_geometries;


    }; 

}

#include "bempp_grid_data_container_impl.hpp"

#endif
