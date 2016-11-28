#include "../../common/common.hpp"
#include "../../common/shared_ptr.hpp"
#include "../../common/boost_make_shared_fwd.hpp"
#include "bempp_triangle_grid.hpp"
#include <dune/geometry/type.hh>
#include <dune/grid/common/entity.hh>
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

        auto data = boost::make_shared<DataContainer>();
        data->init(vertexContainer, elementContainer);
        Dune::Entity<0, 2, const TriangleGrid, EntityImp> entity(EntityImp<0, 2, const TriangleGrid>(data, 0, 0));
        for (int i = 0; i < 3; ++i){
            auto node = entity.subEntity<1>(i)->geometry().center(); 
            for (int j = 0; j < 3; ++j) std::cout << node[j] << " ";
            std::cout << std::endl;
        }
        Geometry<2, 3, const TriangleGrid > geom(geometryType, vertices);
        std::cout << "Number of edges: " << data->numberOfEdges(0) << std::endl;

        auto grid = TriangleGrid(data);

        std::cout << "Test the iterator" << std::endl;

        for (auto it = grid.lbegin<1>(0); it != grid.lend<1>(0); ++it)
            std::cout << it->geometry().center() << std::endl;;

        IdSetImp idSet;
        for (auto it = grid.lbegin<0>(0); it != grid.lend<0>(0); ++it)
            std::cout << idSet.subId(*it, 0, 2) << std::endl;

    }



}
