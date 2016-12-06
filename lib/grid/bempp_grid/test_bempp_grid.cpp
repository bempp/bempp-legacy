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
            boost::make_shared<std::vector<std::vector<unsigned int>>>(2);
        auto& elements = *elementContainer;
        elements[0].resize(3);
        elements[0][0] = 0;
        elements[0][1] = 1;
        elements[0][2] = 2;

        elements[1].resize(3);
        elements[1][0] = 2;
        elements[1][1] = 1;
        elements[1][2] = 3;

        Dune::GridFactory<Dune::Grid<2, 3, double, TriangleGridFamily>> factory;
        for (int i = 0; i < vertices.size(); ++i) 
            factory.insertVertex(vertices[i]);

        Dune::GeometryType type;
        type.makeSimplex(2);
        for (int i = 0; i < elements.size(); ++i)
            factory.insertElement(type, elements[i]);

        shared_ptr<Dune::Grid<2, 3, double, TriangleGridFamily>> grid(factory.createGrid());

        std::cout << "Test the iterator" << std::endl;
        auto view = grid->leafGridView();
        const auto& indexSet = view.indexSet();

    }



}
