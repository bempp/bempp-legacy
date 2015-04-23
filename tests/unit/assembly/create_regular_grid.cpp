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

#include "create_regular_grid.hpp"
#include "grid/grid_factory.hpp"
#include "common/eigen_support.hpp"

using namespace Bempp;

Bempp::shared_ptr<Bempp::Grid> createRegularTriangularGrid(
    int nElementsX, int nElementsY, double maxX, double maxY)
{
    Bempp::GridParameters params;
    params.topology = Bempp::GridParameters::TRIANGULAR;

    const int dimGrid = 2;
    Vector<double> lowerLeft(dimGrid);
    Vector<double> upperRight(dimGrid);
    Vector<int> nElements(dimGrid);
    lowerLeft.fill(0);
    upperRight(0) = maxX;
    upperRight(1) = maxY;
    nElements(0) = nElementsX;
    nElements(1) = nElementsY;

    return Bempp::GridFactory::createStructuredGrid(
        params, lowerLeft, upperRight, nElements);
}
