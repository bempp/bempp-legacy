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

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_view.hpp"
#include "grid/index_set.hpp"
#include "grid/vtk_writer.hpp"

using namespace Bempp;
using std::cout;
using std::endl;

/**
  Load a mesh from a Gmsh file and export a VTK file 'vtkhead.vtu' containing a
  data series 'z_coord' which maps each vertex to its z coordinate. This series
  can be subsequently visulatised in ParaView.
  */
int main()
{
    const char MESH_FNAME[] = "head.gmsh";

    // Import the grid
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

    std::auto_ptr<Grid> grid(GridFactory::importGmshGrid(params, std::string(MESH_FNAME),
                             true, // verbose
                             false)); // insertBoundarySegments

    // Create a leaf view
    std::auto_ptr<GridView> leafGridView(grid->leafView());
    const IndexSet& indexSet = leafGridView->indexSet();

    const int nVertices = leafGridView->entityCount(2);
    arma::Mat<double> data(1,nVertices);
    arma::Col<double> elementCenter;

    cout << "Traversing the grid..." << endl;
    std::auto_ptr<EntityIterator<2> > leafVertexIt(leafGridView->entityIterator<2>());
    while (!leafVertexIt->finished()) {
        const Entity<2>& e = leafVertexIt->entity();

        const Geometry& geo = e.geometry();
        geo.center(elementCenter);
        data(0,indexSet.entityIndex(e)) = elementCenter(2); // the z coordinate

        leafVertexIt->next();
    }

    cout << "Gathering data..." << endl;
    std::auto_ptr<VtkWriter> vtkWriter = leafGridView->vtkWriter();
    vtkWriter->addVertexData(data, "z_coord");
    cout << "Exporting data..." << endl;
    std::string outputFname = vtkWriter->write("vtkhead");
    cout << "Data exported to '" << outputFname << "'." << endl;
}


