#include "p1_data_container.hpp"


namespace Bempp {

    P1DataContainer::P1DataContainer(const shared_ptr<P1DataContainer::NodesContainer>& nodes,
                                     const shared_ptr<P1DataContainer::ElementsContainer>& elements) :
        m_nodes({nodes}), m_elements({elements}), m_levels(0) {}


    const P1DataContainer::NodesContainer& P1DataContainer::nodes(int level) const {
        return *(m_nodes[level]);
    }

    const P1DataContainer::ElementsContainer& P1DataContainer::elements(int level) const{
        return *(m_elements[level]);
    }
    
    void P1DataContainer::populateData(int level) {
    }        
}
