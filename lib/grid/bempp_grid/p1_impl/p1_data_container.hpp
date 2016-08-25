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

            P1DataContainer(const shared_ptr<NodesContainer>& nodes, 
                    const shared_ptr<ElementsContainer>& elements);

            const NodesContainer& nodes(int level) const;
            const ElementsContainer& elements(int level) const;

        private:

            void populateData(int level);

            int m_levels;
            std::vector<shared_ptr<NodesContainer>> m_nodes;
            std::vector<shared_ptr<ElementsContainer>> m_elements;
            std::vector<shared_ptr<EdgesContainer>> m_edges;

    }; 

}


#endif
