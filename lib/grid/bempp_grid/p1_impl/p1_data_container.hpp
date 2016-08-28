#ifndef p1_data_container_hpp
#define p1_data_container_hpp

#include "../../../common/common.hpp"
#include "../../../common/shared_ptr.hpp"
#include <dune/common/fvector.hh>
#include <vector>
#include <array>

namespace Bempp {


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

            int numberOfNodes(int level) const;
            int numberOfElements(int level) const;
            int numberOfEdges(int level) const;
            int levels() const;

        private:

            void populateData(int level);

            int m_levels;
            std::vector<shared_ptr<NodesContainer>> m_nodes;
            std::vector<shared_ptr<ElementsContainer>> m_elements;
            std::vector<shared_ptr<EdgesContainer>> m_edges;
            std::vector<std::vector<std::array<std::size_t, 3>>> m_element2Edges;

    }; 

}


#endif
