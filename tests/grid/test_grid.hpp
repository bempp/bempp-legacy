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

#ifndef bempp_test_grid_hpp
#define bempp_test_grid_hpp

#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"

/** Fixture class for Bempp::Grid tests */
class SimpleTriangularGridManager {
public:
    /** Create two identical simple 2D structured Bempp grids composed of 2 * 3 * 4 triangles.
        Store an auto_ptr to the first one as bemppGrid; from the second one,
        extract the pointer to the underlying Dune grid and store it as duneGrid. */
    SimpleTriangularGridManager() {
        // Create a Bempp grid
        bemppGrid = createGrid();

        // Create an identical Dune grid
        dummyBemppGrid = createGrid();
        Bempp::DefaultGrid* castedBemppGrid = dynamic_cast<Bempp::DefaultGrid*>(&*dummyBemppGrid);
        duneGrid = &castedBemppGrid->duneGrid();
    }

    // No destructor is needed since auto_ptrs release memory automatically

private:
    std::auto_ptr<Bempp::Grid> createGrid() {
        Bempp::GridParameters params;
        params.topology = Bempp::GridParameters::TRIANGULAR;

        const int dimGrid = 2;
        arma::Col<Bempp::ctype> lowerLeft(dimGrid);
        arma::Col<Bempp::ctype> upperRight(dimGrid);
        arma::Col<unsigned int> nElements(dimGrid);
        lowerLeft.fill(0);
        upperRight.fill(1);
        nElements(0) = 3;
        nElements(1) = 4;

        return Bempp::GridFactory::createStructuredGrid(params, lowerLeft, upperRight, nElements);
    }

public:
    typedef Bempp::DefaultDuneGrid DuneGrid;
    std::auto_ptr<Bempp::Grid> bemppGrid;
    DuneGrid* duneGrid;
private:
    std::auto_ptr<Bempp::Grid> dummyBemppGrid; // this grid is used only to manage the life of the duneGrid pointer
};

#endif
