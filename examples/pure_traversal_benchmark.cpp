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

#include <iostream>
#include <memory> // auto_ptr
#include <sys/time.h>

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
    const int dimWorld = 3;
    typedef DefaultGrid::DuneGridType DuneGrid;
    typedef Dune::StructuredGridFactory<DuneGrid> StructGridFactory;
    typedef DuneGrid::ctype ctype;
    const int dimGrid = 2;

    // Construct a Dune structured grid.
    // Currently this procedure is not encapsulated in Bempp.
    // Such encapsulation would be much easier if
    // StructuredGridFactory::createSimplexGrid() returned a standard pointer rather
    // than a shared pointer, because Bempp::ThreeD::Grid doesn't have a place to
    // store the shared pointer.

    const Dune::FieldVector<ctype, dimWorld> lowerLeft(0);
    const Dune::FieldVector<ctype, dimWorld> upperRight(1);
    Dune::array<unsigned int, dimGrid> nElements;
    const int N_ELEMENTS = 1000;
    const int N_TRIALS = 20;

    nElements[0] = N_ELEMENTS;
    nElements[1] = N_ELEMENTS + 1;

    Dune::shared_ptr<DuneGrid> duneGrid =
        StructGridFactory::createSimplexGrid(lowerLeft, upperRight, nElements);

    std::cout << nElements[0] * nElements[1] << " elements created" << std::endl;

    //////////////////// BEMPP OBJECTS ///////////////////

    std::cout << "Using Bempp objects..." << std::endl;

    // Wrap the grid in a Bempp object
    DefaultGrid grid(duneGrid.get());

    // Create a leaf view
    std::auto_ptr<GridView> leafGridView(grid.leafView());

    {
        std::cout << "Iterating over faces..." << std::endl;

        volatile int j = 0;
        timeval start, end;
        gettimeofday(&start, 0);

        for (int i = 0; i < N_TRIALS; ++i) {
            std::auto_ptr<EntityIterator<0> > leafFaceIt = leafGridView->entityIterator<0>();
            while (!leafFaceIt->finished()) {
                ++j;
                leafFaceIt->next();
            }
        }
        gettimeofday(&end, 0);

        double total_time = (end.tv_sec + end.tv_usec / 1000000.) -
                            (start.tv_sec + start.tv_usec / 1000000.);
        std::cout << "Traversed faces: " << j << '\n';
        std::cout << "Traversal time: " << total_time << std::endl;
    }

    {
        std::cout << "Iterating over vertices..." << std::endl;

        volatile int j = 0;
        timeval start, end;
        gettimeofday(&start, 0);

        for (int i = 0; i < N_TRIALS; ++i) {
            std::auto_ptr<EntityIterator<2> > leafVertexIt = leafGridView->entityIterator<2>();
            while (!leafVertexIt->finished()) {
                ++j;
                leafVertexIt->next();
            }
        }
        gettimeofday(&end, 0);

        double total_time = (end.tv_sec + end.tv_usec / 1000000.) -
                            (start.tv_sec + start.tv_usec / 1000000.);
        std::cout << "Traversed vertices: " << j << '\n';
        std::cout << "Traversal time: " << total_time << std::endl;
    }

    //////////////////// DUNE OBJECTS ///////////////////

    std::cout << "Using Dune objects..." << std::endl;

    DefaultDuneGrid::LeafGridView duneLeafGridView = duneGrid->leafView();

    {
        std::cout << "Iterating over faces..." << std::endl;

        volatile int j = 0;
        timeval start, end;
        gettimeofday(&start, 0);

        for (int i = 0; i < N_TRIALS; ++i) {
            typedef DefaultDuneGrid::LeafGridView::Codim<0> Codim;
            Codim::Iterator leafFaceIt = duneLeafGridView.begin<0>();
            Codim::Iterator leafEnd = duneLeafGridView.end<0>();
            for (; leafFaceIt != leafEnd; ++leafFaceIt) {
                ++j;
            }
        }
        gettimeofday(&end, 0);

        double total_time = (end.tv_sec + end.tv_usec / 1000000.) -
                            (start.tv_sec + start.tv_usec / 1000000.);
        std::cout << "Traversed faces: " << j << '\n';
        std::cout << "Traversal time: " << total_time << std::endl;
    }


    {
        std::cout << "Iterating over vertices..." << std::endl;

        volatile int j = 0;
        timeval start, end;
        gettimeofday(&start, 0);

        for (int i = 0; i < N_TRIALS; ++i) {
            typedef DefaultDuneGrid::LeafGridView::Codim<2> Codim;
            Codim::Iterator leafFaceIt = duneLeafGridView.begin<2>();
            Codim::Iterator leafEnd = duneLeafGridView.end<2>();
            for (; leafFaceIt != leafEnd; ++leafFaceIt) {
                ++j;
            }
        }
        gettimeofday(&end, 0);

        double total_time = (end.tv_sec + end.tv_usec / 1000000.) -
                            (start.tv_sec + start.tv_usec / 1000000.);
        std::cout << "Traversed vertices: " << j << '\n';
        std::cout << "Traversal time: " << total_time << std::endl;
    }

}
