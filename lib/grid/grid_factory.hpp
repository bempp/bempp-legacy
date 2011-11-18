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

#ifndef bempp_grid_factory_hpp
#define bempp_grid_factory_hpp

#include "structured_grid_factory.hpp"
#include "common.hpp"

#include <armadillo>
#include <memory>

namespace Bempp {

class Grid;

struct GridParameters
{
    enum Topology {
        LINEAR, /**< one-dimensional grid */
        TRIANGULAR, QUADRILATERAL, HYBRID
    } topology;
};

class GridFactory
{
public:
    static std::auto_ptr<Grid> createStructuredGrid(const GridParameters& params,
                                                    const arma::Col<ctype>& lowerLeft,
                                                    const arma::Col<ctype>& upperRight,
                                                    const arma::Col<unsigned int>& nElements);

    static std::auto_ptr<Grid> importGmshGrid(const GridParameters& params,
                                              const std::string& fileName,
                                              bool verbose=true, bool insertBoundarySegments=true);

    static std::auto_ptr<Grid> importGmshGrid(const GridParameters& params,
                                              const std::string& fileName,
                                              std::vector<int>& boundaryId2PhysicalEntity,
                                              std::vector<int>& elementIndex2PhysicalEntity,
                                              bool verbose=true, bool insertBoundarySegments=true);
};

} // namespace Bempp

#endif
