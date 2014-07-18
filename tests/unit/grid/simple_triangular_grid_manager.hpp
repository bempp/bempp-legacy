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

#ifndef bempp_simple_triangular_grid_manager_hpp
#define bempp_simple_triangular_grid_manager_hpp

#include "grid/dune.hpp"
#include "grid/grid.hpp"
#include "common/shared_ptr.hpp"

/** Fixture class for Bempp::Grid tests */
class SimpleTriangularGridManager
{
public:
    enum { N_ELEMENTS_X = 3, N_ELEMENTS_Y = 4 };
    typedef Bempp::Default2dIn3dDuneGrid DuneGrid;

    /** Create two identical simple 2D structured Bempp grids composed of 2 * 3 * 4 triangles.
        Store an auto_ptr to the first one as bemppGrid; from the second one,
        extract the pointer to the underlying Dune grid and store it as duneGrid. */
    SimpleTriangularGridManager() {
        // Create a Bempp grid
        bemppGrid = createGrid();

        // Create an identical Dune grid
        duneGrid = createDuneGrid();
    }

    // No destructor is needed since auto_ptrs release memory automatically

private:
    Bempp::shared_ptr<Bempp::Grid> createGrid();
    Bempp::shared_ptr<DuneGrid> createDuneGrid();

public:
    Bempp::shared_ptr<Bempp::Grid> bemppGrid;
    Bempp::shared_ptr<DuneGrid> duneGrid;
};

#endif
