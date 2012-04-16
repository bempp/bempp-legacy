// Copyright (C) 2011 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include <armadillo>
#include <iostream>
#include <memory> // auto_ptr

#include "assembly/assembly_options.hpp"
#include "fiber/ordinary_function.hpp"
#include "assembly/grid_function.hpp"

#include "fiber/standard_local_assembler_factory_for_operators_on_surfaces.hpp"

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry_factory.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_view.hpp"
#include "grid/index_set.hpp"
#include "grid/mapper.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include <cmath>

using namespace Bempp;
using std::cout;
using std::endl;

enum MeshVariant
{
    TWO_DISJOINT_TRIANGLES,
    TWO_TRIANGLES_SHARING_VERTEX_0,
    TWO_TRIANGLES_SHARING_VERTICES_2_AND_0,
    TWO_TRIANGLES_SHARING_VERTICES_1_AND_0,
    TWO_TRIANGLES_SHARING_EDGES_0_AND_0,
    TWO_TRIANGLES_SHARING_EDGES_1_AND_0,
    SIMPLE_MESH_9,
    CUBE_12,
    CUBE_384
};

inline bool approxEqual(double x, double y)
{
    return fabs(x - y) / fabs((x + y) / 2.) < 1e-3;
}

class MyFunctor
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef double ValueType;

    // Number of components of the function's argument
    static const int argumentDimension = 3;

    // Number of components of the function's result
    static const int resultDimension = 1;

    // Evaluate the function at the point "point" and store result in
    // the array "result"
    inline void evaluate(const arma::Col<ValueType>& point,
                  arma::Col<ValueType>& result) const {
        // result(0) = sin(0.5 * point(0) * cos(0.25 * point(2))) * cos(point(1));
        result(0) = 1.;
    }
};


/**
    A script for rudimentary testing of the single-layer-potential operator.
  */
int main()
{
    const MeshVariant meshVariant = CUBE_12;

    const char TWO_DISJOINT_TRIANGLES_FNAME[] = "two_disjoint_triangles.msh";
    const char TWO_TRIANGLES_SHARING_VERTEX_0_FNAME[] = "two_triangles_sharing_vertex_0.msh";
    const char TWO_TRIANGLES_SHARING_VERTICES_2_AND_0_FNAME[] = "two_triangles_sharing_vertices_2_and_0.msh";
    const char TWO_TRIANGLES_SHARING_VERTICES_1_AND_0_FNAME[] = "two_triangles_sharing_vertices_1_and_0.msh";
    const char TWO_TRIANGLES_SHARING_EDGES_0_AND_0_FNAME[] = "two_triangles_sharing_edges_0_and_0.msh";
    const char TWO_TRIANGLES_SHARING_EDGES_1_AND_0_FNAME[] = "two_triangles_sharing_edges_1_and_0.msh";
    const char SIMPLE_MESH_9_FNAME[] = "simple_mesh_9_elements.msh";
    const char CUBE_12_FNAME[] = "cube-12.msh";
    const char CUBE_384_FNAME[] = "cube-384.msh";

    const char* MESH_FNAME = 0;
    switch (meshVariant) {
    case TWO_DISJOINT_TRIANGLES:
        MESH_FNAME = TWO_DISJOINT_TRIANGLES_FNAME; break;
    case TWO_TRIANGLES_SHARING_VERTEX_0:
        MESH_FNAME = TWO_TRIANGLES_SHARING_VERTEX_0_FNAME; break;
    case TWO_TRIANGLES_SHARING_VERTICES_2_AND_0:
        MESH_FNAME = TWO_TRIANGLES_SHARING_VERTICES_2_AND_0_FNAME; break;
    case TWO_TRIANGLES_SHARING_VERTICES_1_AND_0:
        MESH_FNAME = TWO_TRIANGLES_SHARING_VERTICES_1_AND_0_FNAME; break;
    case TWO_TRIANGLES_SHARING_EDGES_0_AND_0:
        MESH_FNAME = TWO_TRIANGLES_SHARING_EDGES_0_AND_0_FNAME; break;
    case TWO_TRIANGLES_SHARING_EDGES_1_AND_0:
        MESH_FNAME = TWO_TRIANGLES_SHARING_EDGES_1_AND_0_FNAME; break;
    case SIMPLE_MESH_9:
        MESH_FNAME = SIMPLE_MESH_9_FNAME; break;
    case CUBE_12:
        MESH_FNAME = CUBE_12_FNAME; break;
    case CUBE_384:
        MESH_FNAME = CUBE_384_FNAME; break;
    default:
        throw std::runtime_error("Invalid mesh name");
    }

    // Import the grid
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

    std::auto_ptr<Grid> grid(GridFactory::importGmshGrid(params, std::string(MESH_FNAME),
                             true, // verbose
                             false)); // insertBoundarySegments

    std::cout << "Elements:\n";
    std::auto_ptr<GridView> view = grid->leafView();
    const Mapper& elementMapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& entity = it->entity();
        arma::Mat<double> corners;
        entity.geometry().getCorners(corners);
        std::cout << "Element #" << elementMapper.entityIndex(entity) << ":\n";
        std::cout << corners << '\n';
        it->next();
    }
    std::cout.flush();

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
    assemblyOptions.switchToDense();

    AcaOptions acaOptions;
    acaOptions.eps = 1e-4;
    acaOptions.maximumRank = 10000;
    acaOptions.minimumBlockSize = 2;
    acaOptions.eta = 0.8;
    acaOptions.recompress = true;
//    assemblyOptions.switchToAca(acaOptions);

//    assemblyOptions.switchToSparse();

    assemblyOptions.switchToTbb();
    assemblyOptions.setSingularIntegralCaching(AssemblyOptions::NO);

    Fiber::AccuracyOptions accuracyOptions; // default
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(accuracyOptions);

    MyFunctor functor;
    Fiber::OrdinaryFunction<MyFunctor> function(functor);

    GridFunction<double> gridFunction(space, function, factory, assemblyOptions);
    arma::Col<double> coefficients = gridFunction.coefficients().asArmadilloVector();
    std::cout << "\nGenerated vector of coefficients:\n" << coefficients << std::endl;

    gridFunction.exportToVtk(VtkWriter::VERTEX_DATA,
                             "V", "rhs_assembly");
}



