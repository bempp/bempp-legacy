#include "../../common/common.hpp"
#include "../../common/shared_ptr.hpp"
#include "../../common/boost_make_shared_fwd.hpp"
#include "p1_grid.hpp"
#include <dune/geometry/type.hh>
#include <iostream>

namespace Bempp {

    void test_grid() {

        using namespace BemppGrid;

        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        auto vertexContainer = 
            boost::make_shared<std::vector<Dune::FieldVector<double, 3>>>(4); 
        auto& vertices = *vertexContainer;
        vertices[0][0] = 0;
        vertices[0][1] = 0;
        vertices[0][2] = 0;

        vertices[1][0] = 1;
        vertices[1][1] = 0;
        vertices[1][2] = 0;

        vertices[2][0] = 0;
        vertices[2][1] = 1;
        vertices[2][2] = 0;

        vertices[3][0] = 1;
        vertices[3][1] = 1;
        vertices[3][2] = 0;

        auto elementContainer = 
            boost::make_shared<std::vector<std::array<std::size_t, 3>>>(2);
        auto& elements = *elementContainer;
        elements[0][0] = 0;
        elements[0][1] = 1;
        elements[0][2] = 2;

        elements[1][0] = 2;
        elements[1][1] = 1;
        elements[1][2] = 3;

        auto data = boost::make_shared<P1DataContainer>();
        data->addLevel(vertexContainer, elementContainer);
        P1EntityImp<0, 2, P1Grid> entity(data, 0, 0); 
        P1GridGeometry<2, 3, P1Grid > geom(geometryType, vertices);
        std::cout << "Number of edges: " << data->numberOfEdges(0) << std::endl;


    }



}
