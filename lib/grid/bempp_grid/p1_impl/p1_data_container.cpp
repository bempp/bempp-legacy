#include "p1_data_container.hpp"
#include <cassert>
#include <utility>

namespace BemppGrid {

    P1DataContainer::P1DataContainer() :
        m_levels(0), m_idCounter(0) {};

    void P1DataContainer::addLevel(const shared_ptr<P1DataContainer::NodesContainer>& nodes,
                                     const shared_ptr<P1DataContainer::ElementsContainer>& elements){

        int nelements = elements->size();
        int nnodes = nodes->size();

        m_nodes.push_back(nodes);
        m_elements.push_back(elements);
        m_edges.push_back(EdgesContainer());
        auto& edges = m_edges[m_levels];
        m_element2Edges.push_back(std::vector<std::array<std::size_t, 3>>());
        auto& element2Edges = m_element2Edges[m_levels];

        element2Edges.resize(nelements);
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> nodes2EdgeIndexPair;
        nodes2EdgeIndexPair.resize(nodes->size());
        for (int elementIndex = 0; elementIndex < nelements; ++elementIndex) {
            const auto& element = (*elements)[elementIndex];
            for (int i = 3; i > 0; --i) { 
                // Strange counting due to edge numbering in Dune
                std::size_t n0 = element[i-1];
                std::size_t n1 = element[i % 3];
                if (n1 > n0) std::swap(n0, n1); // Number edges from smaller to larger vertex
                // Check if vertex already exists
                bool edgeExists = false;
                std::size_t edgeIndex;
                for (const auto& indexPair: nodes2EdgeIndexPair[n0])
                    if (n1 == indexPair.first) {
                        edgeExists = true;
                        edgeIndex = indexPair.second;
                        break;
                    }
                if (!edgeExists){
                    // Create Edge
                    edgeIndex = edges.size();
                    nodes2EdgeIndexPair[n0].push_back(std::pair<std::size_t, std::size_t>(n1, edgeIndex));
                    edges.push_back(std::array<std::size_t, 2>({n0, n1}));
                }
                element2Edges[elementIndex][i-1] = edgeIndex;
            }
        }

        // Fill the connectivity data arrays
        m_edge2Elements.push_back(std::vector<std::vector<std::size_t>>());
        m_node2Elements.push_back(std::vector<std::vector<std::size_t>>());
        m_node2Edges.push_back(std::vector<std::vector<std::size_t>>());
        
        auto& edge2Elements = m_edge2Elements[m_levels];
        auto& node2Elements = m_node2Elements[m_levels];
        auto& node2Edges = m_node2Edges[m_levels];

        edge2Elements.resize(edges.size());
        node2Elements.resize(nnodes);
        node2Edges.resize(nnodes);

        for (std::size_t i = 0; i < nelements; ++i){
            for (int j = 0; j < 3; ++j){
                edge2Elements[element2Edges[i][j]].push_back(i);
            }
        }

        for (std::size_t i = 0; i < nelements; ++i)
            for (int j = 0; j < 3; ++j)
                node2Elements[(*elements)[i][j]].push_back(i);

        for (std::size_t i = 0; i < edges.size(); ++i)
            for (int j = 0; j < 2; ++j)
                node2Edges[edges[i][j]].push_back(i);

        // Now generate the Ids

        m_nodeIds.push_back(std::vector<std::size_t>(m_nodes[m_levels]->size()));
        m_elementIds.push_back(std::vector<std::size_t>(m_elements[m_levels]->size()));
        m_edgeIds.push_back(std::vector<std::size_t>(m_edges[m_levels].size()));

        for (std::size_t i = 0; i < m_nodes[m_levels]->size(); ++i)
            m_nodeIds[m_levels][i] = m_idCounter++;

        for (std::size_t i = 0; i < m_elements[m_levels]->size(); ++i)
            m_elementIds[m_levels][i] = m_idCounter++;

        for (std::size_t i = 0; i < m_edges[m_levels].size(); ++i)
            m_edgeIds[m_levels][i] = m_idCounter++;

        m_levels++;
    }


    const P1DataContainer::NodesContainer& P1DataContainer::nodes(int level) const {
        assert(level < m_levels);
        return *(m_nodes[level]);
    }

    const P1DataContainer::ElementsContainer& P1DataContainer::elements(int level) const{
        assert(level < m_levels);
        return *(m_elements[level]);
    }

    const P1DataContainer::EdgesContainer& P1DataContainer::edges(int level) const {
        assert(level < m_levels);
        return m_edges[level];
    }

    int P1DataContainer::numberOfNodes(int level) const {

        assert(level < m_levels);
        return m_nodes[level]->size();

    }
   
    int P1DataContainer::numberOfElements(int level) const {

        assert(level < m_levels);
        return m_elements[level]->size();

    }

    int P1DataContainer::numberOfEdges(int level) const {

        assert(level < m_levels);
        return m_edges[level].size();

    }


    const std::array<std::size_t, 3>& P1DataContainer::element2Edges(
            int level, std::size_t elementIndex) const {
        assert(level < m_levels);
        return m_element2Edges[level][elementIndex];

    }

    const std::vector<size_t>& P1DataContainer::edge2Elements(int level, std::size_t edgeIndex) const {
        assert(level < m_levels);
        return m_edge2Elements[level][edgeIndex];

    }
    const std::vector<size_t>& P1DataContainer::node2Elements(int level, std::size_t nodeIndex) const {
        assert(level < m_levels);
        return m_node2Elements[level][nodeIndex];

    }

    const std::vector<size_t>& P1DataContainer::node2Edges(int level, std::size_t nodeIndex) const {
        assert(level < m_levels);
        return m_node2Edges[level][nodeIndex];

    }

}
