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

#include "grid/entity_iterator.hpp"
#include "grid/entity.hpp"
#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid_view.hpp"
#include "grid/geometry.hpp"
#include "grid/mapper.hpp"

using namespace Bempp;

shared_ptr<Grid>
loadTriangularMeshFromFile(const char* fileName)
{
    // Import the grid
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;

    return GridFactory::importGmshGrid(params, std::string(fileName));
}

void dumpElementList(const Grid* grid)
{
    std::cout << "Elements:\n";
    std::unique_ptr<GridView> view = grid->leafView();
    const Mapper& elementMapper = view->elementMapper();
    std::unique_ptr<EntityIterator<0> > it = view->entityIterator<0>();
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
