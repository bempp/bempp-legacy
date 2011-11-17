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

/* NOTE: the version included with Dune 2.1 has a bug. The declaration of
StructuredGridFactory::insertVertices() needs to be changed to

static void insertVertices(GridFactory<GridType>& factory,
                           const FieldVector<ctype,dimworld>& lowerLeft,
                           const FieldVector<ctype,dimworld>& upperRight,
                           const array<unsigned int,dim>& vertices)

(note the change of dim to dimworld in two places).
*/
#include <dune/grid/utility/structuredgridfactory.hh>

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"

using namespace Bempp;

int main()
{
    // Grid type and dimensions
    typedef DefaultGrid::DuneGridType DuneGrid;
    typedef Dune::StructuredGridFactory<DuneGrid> StructGridFactory;
    typedef DuneGrid::ctype ctype;
    const int dimGrid = 2;
    const int dimWorld = 3;

    // Construct a Dune structured grid.
    // Currently this procedure is not encapsulated in Bempp.
    // Such encapsulation would be much easier if
    // StructuredGridFactory::createSimplexGrid() returned a standard pointer rather
    // than a shared pointer, because Bempp::ThreeD::Grid doesn't have a place to
    // store the shared pointer.

    const Dune::FieldVector<ctype, dimWorld> lowerLeft(0);
    const Dune::FieldVector<ctype, dimWorld> upperRight(1);
    Dune::array<unsigned int, dimGrid> nElements;
    nElements[0] = 2;
    nElements[1] = 3;

    Dune::shared_ptr<DuneGrid> duneGrid =
        StructGridFactory::createSimplexGrid(lowerLeft, upperRight, nElements);

    // Wrap the grid in a Bempp object
    DefaultGrid grid(duneGrid.get());

    // Create a leaf view
    std::auto_ptr<GridView> leafGridView(grid.leafView());

    {
        std::cout << "Iterating over faces..." << std::endl;
        std::auto_ptr<EntityIterator<0> > leafFaceIt(leafGridView->entityIterator<0>());
        while (!leafFaceIt->finished()) {
            const Entity<0>& e = leafFaceIt->entity();
            Dune::GeometryType gt = e.type();
            const Geometry& geo = e.geometry();
            // Dune::FieldVector<ctype, dimWorld> elementCenter;
            arma::Col<double> elementCenter;
            geo.center(elementCenter);

            std::cout << "  visiting leaf element " << gt
                      << " with centre at " << elementCenter << std::endl;

            leafFaceIt->next();
        }
    }

    {
        std::cout << "Iterating over vertices..." << std::endl;
        std::auto_ptr<EntityIterator<2> > leafVertexIt(leafGridView->entityIterator<2>());
        while (!leafVertexIt->finished()) {
            const Entity<2>& e = leafVertexIt->entity();
            Dune::GeometryType gt = e.type();
            const Geometry& geo = e.geometry();
            // Dune::FieldVector<ctype, dimWorld> elementCenter = geo.center();
            arma::Col<double> elementCenter;
            geo.center(elementCenter);

            std::cout << "  visiting leaf element " << gt
                      << " with centre at " << elementCenter << std::endl;

            leafVertexIt->next();
        }
    }

    {
        std::cout << "Iterating over vertices of the second face" << std::endl;
        std::auto_ptr<EntityIterator<0> > leafFaceIt(leafGridView->entityIterator<0>());
        leafFaceIt->next();
        std::auto_ptr<EntityIterator<2> > vertexIt(leafFaceIt->entity().subEntityIterator<2>());
        while (!vertexIt->finished()) {
            const Entity<2>& e = vertexIt->entity();
            Dune::GeometryType gt = e.type();
            const Geometry& geo = e.geometry();
            arma::Col<double> elementCenter;
            geo.center(elementCenter);

            std::cout << "  visiting vertex " << gt
                      << " with centre at " << elementCenter << std::endl;

            vertexIt->next();
        }
    }

    {
        std::cout << "Now the program will throw an exception!" << std::endl; 
        std::auto_ptr<EntityIterator<3> > invalidIt(leafGridView->entityIterator<3>());
    }
}
