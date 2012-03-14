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
#include "assembly/identity_operator.hpp"
#include "assembly/single_layer_potential_3d.hpp"

#include "assembly/discrete_scalar_valued_linear_operator.hpp"
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
#include "grid/vtk_writer.hpp"

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

enum OperatorVariant
{
    SINGLE_LAYER_POTENTIAL,
    IDENTITY
};

inline bool approxEqual(double x, double y)
{
    return fabs(x - y) / fabs((x + y) / 2.) < 1e-3;
}

/**
    A script for rudimentary testing of the single-layer-potential operator.
  */
int main()
{
    const MeshVariant meshVariant = TWO_TRIANGLES_SHARING_VERTEX_0;
    const OperatorVariant opVariant = SINGLE_LAYER_POTENTIAL;

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
        entity.geometry().corners(corners);
        std::cout << "Element #" << elementMapper.entityIndex(entity) << ":\n";
        std::cout << corners << '\n';
        it->next();
    }
    std::cout.flush();

    PiecewiseLinearContinuousScalarSpace<double> space(*grid);
    space.assignDofs();

    AssemblyOptions assemblyOptions;
    assemblyOptions.mode = ASSEMBLY_MODE_DENSE;
    // assemblyOptions.mode = ASSEMBLY_MODE_ACA;
    // assemblyOptions.mode = ASSEMBLY_MODE_SPARSE;
    assemblyOptions.acaEps = 1e-4;
    assemblyOptions.acaMaximumRank = 10000;
    assemblyOptions.acaMinimumBlockSize = 2;
    assemblyOptions.acaEta = 0.8;
    assemblyOptions.acaRecompress = true;

    Fiber::OpenClOptions openClOptions;
    openClOptions.useOpenCl = false;
    Fiber::StandardLocalAssemblerFactoryForOperatorsOnSurfaces<double, GeometryFactory>
            factory(openClOptions);

    typedef std::auto_ptr<LinearOperator<double> > LinearOperatorPtr;
    LinearOperatorPtr op;
    switch (opVariant)
    {
    case SINGLE_LAYER_POTENTIAL:
        op = LinearOperatorPtr(new SingleLayerPotential3D<double>); break;
    case IDENTITY:
        op = LinearOperatorPtr(new IdentityOperator<double>); break;
    default:
        throw std::runtime_error("Invalid operator");
    }

    std::auto_ptr<DiscreteScalarValuedLinearOperator<double> > result =
            op->assembleWeakForm(space, space, factory, assemblyOptions);

    arma::Mat<double> resultMatrix = result->asMatrix();
    std::cout << "\nGenerated matrix:\n" << resultMatrix << std::endl;

    if (opVariant == SINGLE_LAYER_POTENTIAL)
    {
        if (meshVariant == TWO_DISJOINT_TRIANGLES)
        {
            // disjoint elements
            assert(approxEqual(resultMatrix(0, 3), 0.0101397));
            assert(approxEqual(resultMatrix(0, 4), 0.00852427));
            assert(approxEqual(resultMatrix(1, 3), 0.0127044));
            // TODO: calculate with reasonable accuracy the result for coincident triangles
            // and make a test
        }
        else if (meshVariant == TWO_TRIANGLES_SHARING_VERTEX_0)
        {
            // elements sharing vertex
            assert(approxEqual(resultMatrix(1, 3), 0.0109119));
            assert(approxEqual(resultMatrix(2, 4), 0.0123149));
        }
        else if (meshVariant == TWO_TRIANGLES_SHARING_VERTICES_2_AND_0)
        {
            // elements sharing vertex
            assert(approxEqual(resultMatrix(0, 3), 0.0109119));
            assert(approxEqual(resultMatrix(1, 4), 0.0123149));
        }
        else if (meshVariant == TWO_TRIANGLES_SHARING_VERTICES_1_AND_0)
        {
            // elements sharing vertex
            assert(approxEqual(resultMatrix(0, 3), 0.0109119));
            assert(approxEqual(resultMatrix(2, 4), 0.0123149));
        }
        else if (meshVariant == TWO_TRIANGLES_SHARING_EDGES_0_AND_0)
        {
            // elements sharing edge
            assert(approxEqual(resultMatrix(2, 3), 0.00861375));
        }
        else if (meshVariant == TWO_TRIANGLES_SHARING_EDGES_1_AND_0)
        {
            // elements sharing edge
            assert(approxEqual(resultMatrix(0, 3), 0.00861375));
        }
    }
    else if (opVariant == IDENTITY)
    {
        if (meshVariant == TWO_DISJOINT_TRIANGLES)
        {
            // disjoint elements
            assert(approxEqual(resultMatrix(0, 0), 0.25));
            assert(approxEqual(resultMatrix(0, 1), 0.125));
            assert(approxEqual(resultMatrix(0, 2), 0.125));
            assert(approxEqual(resultMatrix(4, 4), 0.5));
            assert(approxEqual(resultMatrix(4, 5), 0.25));
        }
        else if (meshVariant == TWO_TRIANGLES_SHARING_VERTEX_0)
        {
            // elements sharing vertex
            assert(approxEqual(resultMatrix(0, 0), 0.75));
            assert(approxEqual(resultMatrix(0, 1), 0.125));
            assert(approxEqual(resultMatrix(0, 2), 0.125));
            assert(approxEqual(resultMatrix(0, 3), 0.25));
            assert(approxEqual(resultMatrix(0, 4), 0.25));
        }
    }
}



