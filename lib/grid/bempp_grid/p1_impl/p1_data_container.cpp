#include "p1_data_container.hpp"
#include <cassert>
#include <utility>

namespace Bempp {

    P1DataContainer::P1DataContainer() :
        m_levels(0) {};

    void P1DataContainer::addLevel(const shared_ptr<P1DataContainer::NodesContainer>& nodes,
                                     const shared_ptr<P1DataContainer::ElementsContainer>& elements){

        m_nodes.push_back(nodes);
        m_elements.push_back(elements);
        int numberOfElements = elements->size();
        auto edges = shared_ptr<EdgesContainer>(
                    new EdgesContainer());

        std::vector<std::array<std::size_t, 3>> element2Edges;
        element2Edges.resize(elements->size());
        std::vector<std::vector<std::pair<std::size_t, std::size_t>>> nodes2EdgeIndexPair;
        nodes2EdgeIndexPair.resize(nodes->size());
        for (int elementIndex = 0; elementIndex < numberOfElements; ++elementIndex) {
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
                    edgeIndex = edges->size();
                    nodes2EdgeIndexPair[n0].push_back(std::pair<std::size_t, std::size_t>(n1, edgeIndex));
                    edges->push_back(std::array<std::size_t, 2>({n0, n1}));
                }
                element2Edges[elementIndex][i-1] = edgeIndex;
            }
        }

        m_edges.push_back(edges);
        m_levels += 1;

    }
                    




    const P1DataContainer::NodesContainer& P1DataContainer::nodes(int level) const {
        assert(level < m_levels);
        return *(m_nodes[level]);
    }

    const P1DataContainer::ElementsContainer& P1DataContainer::elements(int level) const{
        assert(level < m_levels);
        return *(m_elements[level]);
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
        return m_edges[level]->size();

    }


}
