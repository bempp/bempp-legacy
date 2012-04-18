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

#include "meshes.hpp"

#include <iostream>

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/mapper.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/entity.hpp"
#include "grid/grid_view.hpp"
#include "grid/index_set.hpp"
#include "grid/vtk_writer.hpp"
#include "grid/geometry.hpp"

using namespace Bempp;

std::auto_ptr<Grid>
loadMesh(MeshVariant mv)
{
    const char TWO_DISJOINT_TRIANGLES_FNAME[] = "two_disjoint_triangles.msh";
    const char TWO_TRIANGLES_SHARING_VERTEX_0_FNAME[] = "two_triangles_sharing_vertex_0.msh";
    const char TWO_TRIANGLES_SHARING_VERTICES_2_AND_0_FNAME[] = "two_triangles_sharing_vertices_2_and_0.msh";
    const char TWO_TRIANGLES_SHARING_VERTICES_1_AND_0_FNAME[] = "two_triangles_sharing_vertices_1_and_0.msh";
    const char TWO_TRIANGLES_SHARING_EDGES_0_AND_0_FNAME[] = "two_triangles_sharing_edges_0_and_0.msh";
    const char TWO_TRIANGLES_SHARING_EDGES_1_AND_0_FNAME[] = "two_triangles_sharing_edges_1_and_0.msh";
    const char SIMPLE_MESH_9_FNAME[] = "simple_mesh_9_elements.msh";
    const char CUBE_12_FNAME[] = "cube-12.msh";
    const char CUBE_12_REORIENTED_FNAME[] = "cube-12-reoriented.msh";
    const char CUBE_384_FNAME[] = "cube-384.msh";
    const char SPHERE_614_FNAME[] = "sphere-614.msh";
    const char SPHERE_2590_FNAME[] = "sphere-2590.msh";
    const char CUBE_6144_FNAME[] = "cube-6144.msh";
    const char CUBE_24576_FNAME[] = "cube-24576.msh";
    const char SPHERE_152_FNAME[] = "sphere-152.msh";
    const char SPHERE_644_FNAME[] = "sphere-644.msh";
    const char SPHERE_41440_FNAME[] = "sphere-41440.msh";
    const char DOUBLE_SPHERE_5162_FNAME[] = "double-sphere-5162.msh";
    const char DOUBLE_SPHERE_5162_FNAME[] = "double-sphere-20624.msh";
    const char SPHERICAL_SHELL_INNER_SURFACE_FNAME[] =
            "spherical_shell_2194_nodes_inner_surface.msh";

    const char* MESH_FNAME = 0;
    switch (mv) {
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
    case CUBE_12_REORIENTED:
        MESH_FNAME = CUBE_12_REORIENTED_FNAME; break;
    case CUBE_384:
        MESH_FNAME = CUBE_384_FNAME; break;
    case SPHERE_614:
    	MESH_FNAME = SPHERE_614_FNAME; break;
    case SPHERE_2590:
    	MESH_FNAME = SPHERE_2590_FNAME; break;
    case CUBE_6144:
        MESH_FNAME = CUBE_6144_FNAME; break;
    case CUBE_24576:
        MESH_FNAME = CUBE_24576_FNAME; break;
    case SPHERE_152:
        MESH_FNAME = SPHERE_152_FNAME; break;
    case SPHERE_644:
        MESH_FNAME = SPHERE_644_FNAME; break;
    case DOUBLE_SPHERE_5162:
        MESH_FNAME = DOUBLE_SPHERE_5162_FNAME; break;
    case DOUBLE_SPHERE_20624:
        MESH_FNAME = DOUBLE_SPHERE_20624_FNAME; break;
    case SPHERE_41440:
        MESH_FNAME = SPHERE_41440_FNAME; break;
    case SPHERICAL_SHELL_INNER_SURFACE:
        MESH_FNAME = SPHERICAL_SHELL_INNER_SURFACE_FNAME; break;
    default:
        throw std::runtime_error("Invalid mesh name");
    }

    // Import the grid
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

    return GridFactory::importGmshGrid(params, std::string(MESH_FNAME),
                                       true, // verbose
                                       false); // insertBoundarySegments
}

void dumpElementList(const Grid* grid)
{
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
}
