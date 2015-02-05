// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "grid_factory.hpp"
#include "concrete_grid.hpp"
#include "dune.hpp"
#include "structured_grid_factory.hpp"

#include "../common/to_string.hpp"

#include <dune/grid/io/file/gmshreader.hh>
#include <stdexcept>
#include <string>
#include <vector>

namespace Bempp {

// Default grid typedef
typedef ConcreteGrid<Default2dIn3dDuneGrid> Default2dIn3dGrid;

shared_ptr<Grid>
GridFactory::createStructuredGrid(const GridParameters &params,
                                  const arma::Col<double> &lowerLeft,
                                  const arma::Col<double> &upperRight,
                                  const arma::Col<unsigned int> &nElements) {
  // TODO: Support quadrilateral and linear grids

  // Check arguments
  if (params.topology != GridParameters::TRIANGULAR)
    throw std::invalid_argument(
        "GridFactory::createStructuredGrid(): unsupported grid topology");
  const int dimGrid = 2;
  if ((int)lowerLeft.n_rows != dimGrid)
    throw std::invalid_argument("GridFactory::createStructuredGrid(): "
                                "lowerLeft must be two-dimensional");
  if ((int)upperRight.n_rows != dimGrid)
    throw std::invalid_argument("GridFactory::createStructuredGrid(): "
                                "lowerLeft must be two-dimensional");
  if ((int)nElements.n_rows != dimGrid)
    throw std::invalid_argument("GridFactory::createStructuredGrid(): "
                                "nElements must be two-dimensional");
  if (nElements(0) < 1 || nElements(1) < 1)
    throw std::invalid_argument("GridFactory::createStructuredGrid(): all "
                                "entries in nElements must be positive");

  // Convert Armadillo vectors to Dune vectors
  // TODO: Write nice conversion functions
  typedef Default2dIn3dDuneGrid::ctype ctype;
  Dune::FieldVector<ctype, dimGrid> duneLowerLeft;
  duneLowerLeft[0] = lowerLeft(0);
  duneLowerLeft[1] = lowerLeft(1);
  Dune::FieldVector<ctype, dimGrid> duneUpperRight;
  duneUpperRight[0] = upperRight(0);
  duneUpperRight[1] = upperRight(1);
  Dune::array<int, dimGrid> duneNElements;
  duneNElements[0] = nElements(0);
  duneNElements[1] = nElements(1);

  std::unique_ptr<Default2dIn3dDuneGrid> apDuneGrid;
  // TODO: Support quadrilateral grids using createCubeGrid()
  apDuneGrid = Dune::BemppStructuredGridFactory<
      Default2dIn3dDuneGrid>::createSimplexGrid(duneLowerLeft, duneUpperRight,
                                                duneNElements);
  return shared_ptr<Grid>(
      new Default2dIn3dGrid(apDuneGrid.release(), GridParameters::TRIANGULAR,
                            true)); // true -> owns Dune grid
}

shared_ptr<Grid> GridFactory::importGmshGrid(const GridParameters &params,
                                             const std::string &fileName,
                                             bool verbose,
                                             bool insertBoundarySegments) {
  std::vector<int> boundaryId2PhysicalEntity;
  std::vector<int> elementIndex2PhysicalEntity;

  if (params.topology == GridParameters::TRIANGULAR) {
    Default2dIn3dDuneGrid *duneGrid =
        Dune::GmshReader<Default2dIn3dDuneGrid>::read(
            fileName, boundaryId2PhysicalEntity, elementIndex2PhysicalEntity,
            verbose, insertBoundarySegments);
    return shared_ptr<Grid>(new Default2dIn3dGrid(
        duneGrid, params.topology, elementIndex2PhysicalEntity,
        true)); // true -> owns Dune grid
  }
  else
    throw std::invalid_argument("GridFactory::importGmshGrid(): "
                                "unsupported grid topology");
}

shared_ptr<Grid>
GridFactory::importGmshGrid(const GridParameters &params,
                            const std::string &fileName,
                            std::vector<int> &boundaryId2PhysicalEntity,
                            std::vector<int> &elementIndex2PhysicalEntity,
                            bool verbose, bool insertBoundarySegments) {
  if (params.topology == GridParameters::TRIANGULAR) {
    Default2dIn3dDuneGrid *duneGrid =
        Dune::GmshReader<Default2dIn3dDuneGrid>::read(
            fileName, boundaryId2PhysicalEntity, elementIndex2PhysicalEntity,
            verbose, insertBoundarySegments);
    return shared_ptr<Grid>(new Default2dIn3dGrid(
        duneGrid, params.topology, elementIndex2PhysicalEntity,
        true)); // true -> owns Dune grid
  }
  else
    throw std::invalid_argument("GridFactory::importGmshGrid(): "
                                "unsupported grid topology");
}

shared_ptr<Grid> GridFactory::createGridFromConnectivityArrays(
    const GridParameters &params, const arma::Mat<double> &vertices,
    const arma::Mat<int> &elementCorners,
    const std::vector<int> &domainIndices) {
  const int dimGrid = 2, dimWorld = 3;
  if (params.topology != GridParameters::TRIANGULAR)
    throw std::invalid_argument("createGridFromConnectivityArrays(): "
                                "unsupported grid topology");
  if (vertices.n_rows != dimWorld)
    throw std::invalid_argument("createGridFromConnectivityArrays(): "
                                "the 'vertices' array "
                                "must have exactly 3 rows");
  if (elementCorners.n_rows < 3)
    throw std::invalid_argument("createGridFromConnectivityArrays(): "
                                "the 'elementCorners' array "
                                "must have at least 3 rows");
  if (!domainIndices.empty() && domainIndices.size() != elementCorners.n_cols)
    throw std::invalid_argument(
        "createGridFromConnectivityArrays(): "
        "'domainIndices' must either be empty or contain as many "
        "elements as 'elementCorners' has columns");

  Dune::GridFactory<Default2dIn3dDuneGrid> factory;

  for (size_t i = 0; i < vertices.n_cols; ++i) {
    Dune::FieldVector<double, dimWorld> v;
    v[0] = vertices(0, i);
    v[1] = vertices(1, i);
    v[2] = vertices(2, i);
    factory.insertVertex(v);
  }

  const GeometryType type(GeometryType::simplex, dimGrid);
  const size_t vertexCount = vertices.n_cols;
  std::vector<unsigned int> corners(3);
  for (size_t i = 0; i < elementCorners.n_cols; ++i) {
    if (elementCorners(0, i) < 0 || elementCorners(0, i) >= vertexCount ||
        elementCorners(1, i) < 0 || elementCorners(1, i) >= vertexCount ||
        elementCorners(2, i) < 0 || elementCorners(2, i) >= vertexCount)
      throw std::invalid_argument("createGridFromConnectivityArrays(): invalid "
                                  "vertex index in element #" +
                                  toString(i));
    corners[0] = elementCorners(0, i);
    corners[1] = elementCorners(1, i);
    corners[2] = elementCorners(2, i);
    factory.insertElement(type, corners);
  }
  shared_ptr<Grid> result;
  if (domainIndices.empty())
    result.reset(new Default2dIn3dGrid(factory.createGrid(),
                                       GridParameters::TRIANGULAR,
                                       true /* owns Dune grid */));
  else
    result.reset(
        new Default2dIn3dGrid(factory.createGrid(), GridParameters::TRIANGULAR,
                              domainIndices, true /* owns Dune grid */));
  return result;
}

} // namespace Bempp
