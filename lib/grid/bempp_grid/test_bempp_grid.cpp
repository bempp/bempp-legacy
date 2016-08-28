#include "../../common/common.hpp"
#include "../../common/shared_ptr.hpp"
#include "../../common/boost_make_shared_fwd.hpp"
#include "p1_grid.hpp"
#include <dune/geometry/type.hh>
#include <iostream>

namespace Bempp {

    void create_bempp_geometry() {

        Dune::GeometryType geometryType;
        geometryType.makeTriangle();
        auto vertexContainer = 
            boost::make_shared<std::vector<Dune::FieldVector<double, 3>>>(3); 
        auto& vertices = *vertexContainer;
        vertices.resize(3);
        vertices[0][0] = 0;
        vertices[0][1] = 0;
        vertices[0][2] = 0;

        vertices[1][0] = 1;
        vertices[1][1] = 0;
        vertices[1][2] = 0;

        vertices[2][0] = 0;
        vertices[2][1] = 1;
        vertices[2][2] = 0;

        auto elementContainer = 
            boost::make_shared<std::vector<std::array<std::size_t, 3>>>(1);
        auto& elements = *elementContainer;
        elements.resize(1);
        elements[0][0] = 0;
        elements[0][1] = 1;
        elements[0][2] = 2;


        auto data = boost::make_shared<P1DataContainer>();
        data->addLevel(vertexContainer, elementContainer);
        P1Entity<0, 2, P1Grid> entity(data, 0, 0); 
        P1GridGeometry<2, 3, P1Grid > geom(geometryType, vertices);
        std::cout << "Geometry created" << std::endl;

    }



}
