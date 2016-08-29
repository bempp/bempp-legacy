#ifndef p1_data_container_hpp
#define p1_data_container_hpp

#include "../../../common/common.hpp"
#include "../bempp_grid_types.hpp"
#include <dune/common/fvector.hh>
#include <vector>
#include <array>

namespace BemppGrid {

    class P1DataContainer {

        public:

            typedef std::vector<Dune::FieldVector<double, 3>> NodesContainer;
            typedef std::vector<std::array<std::size_t, 3>> ElementsContainer;
            typedef std::vector<std::array<std::size_t, 2>> EdgesContainer;

            P1DataContainer();

            void addLevel( const shared_ptr<NodesContainer>& nodes, 
                    const shared_ptr<ElementsContainer>& elements);

            const NodesContainer& nodes(int level) const;
            const ElementsContainer& elements(int level) const;
            const EdgesContainer& edges(int level) const;

            int numberOfNodes(int level) const;
            int numberOfElements(int level) const;
            int numberOfEdges(int level) const;
            int levels() const;

            const std::vector<size_t>& edge2Elements(int level, std::size_t edgeIndex) const;
            const std::vector<size_t>& node2Elements(int level, std::size_t nodeIndex) const;
            const std::vector<size_t>& node2Edges(int level, std::size_t nodeIndex) const; 

        private:

            void populateData(int level);

            int m_levels;
            std::vector<shared_ptr<NodesContainer>> m_nodes;
            std::vector<shared_ptr<ElementsContainer>> m_elements;
            std::vector<EdgesContainer> m_edges;
            std::vector<std::vector<std::array<std::size_t, 3>>> m_element2Edges;
            std::vector<std::vector<std::vector<std::size_t>>> m_edge2Elements;
            std::vector<std::vector<std::vector<std::size_t>>> m_node2Elements;
            std::vector<std::vector<std::vector<std::size_t>>> m_node2Edges;

    }; 

}


#endif
