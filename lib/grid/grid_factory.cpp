#include "grid_factory.hpp"
#include "grid.hpp"
#include "structured_grid_factory.hpp"

#include <dune/grid/io/file/gmshreader.hh>
#include <stdexcept>
#include <vector>

namespace Bempp {

std::auto_ptr<Grid> GridFactory::createStructuredGrid(
        const GridParameters& params, const arma::Col<ctype>& lowerLeft,
        const arma::Col<ctype>& upperRight, const arma::Col<unsigned int>& nElements)
{
    /** \todo Support quadrilateral and linear grids */

    // Check arguments
    if (params.topology != GridParameters::TRIANGULAR)
        throw std::invalid_argument("GridFactory::createStructuredGrid(): unsupported grid topology");
    const int dimGrid = 2;
    if (lowerLeft.n_rows != dimGrid)
        throw std::invalid_argument("GridFactory::createStructuredGrid(): lowerLeft must be two-dimensional");
    if (upperRight.n_rows != dimGrid)
        throw std::invalid_argument("GridFactory::createStructuredGrid(): lowerLeft must be two-dimensional");
    if (nElements.n_rows != dimGrid)
        throw std::invalid_argument("GridFactory::createStructuredGrid(): nElements must be two-dimensional");
    if (nElements(0) < 1 || nElements(1) < 1)
        throw std::invalid_argument("GridFactory::createStructuredGrid(): all entries in nElements must be positive");

    // Convert Armadillo vectors to Dune vectors
    /** \todo Write nice conversion functions */
    Dune::FieldVector<ctype,dimGrid> duneLowerLeft;
    duneLowerLeft[0] = lowerLeft(0); duneLowerLeft[1] = lowerLeft(1);
    Dune::FieldVector<ctype,dimGrid> duneUpperRight;
    duneUpperRight[0] = upperRight(0); duneUpperRight[1] = upperRight(1);
    Dune::array<unsigned int,dimGrid> duneNElements;
    duneNElements[0] = nElements(0); duneNElements[1] = nElements(1);

    std::auto_ptr<DefaultDuneGrid> apDuneGrid;
    /** \todo Support quadrilateral grids using createCubeGrid() */
    apDuneGrid = Dune::BemppStructuredGridFactory<DefaultDuneGrid>::
            createSimplexGrid(duneLowerLeft, duneUpperRight, duneNElements);
    return std::auto_ptr<Grid>(new DefaultGrid(apDuneGrid.release(), true)); // true -> owns Dune grid
}

std::auto_ptr<Grid> importGmshGrid(const GridParameters& params,
                                   const std::string& fileName,
                                   bool verbose=true, bool insertBoundarySegments=true)
{
    // Check arguments
    if (params.topology != GridParameters::TRIANGULAR)
        throw std::invalid_argument("GridFactory::importGmshGrid(): unsupported grid topology");

    DefaultDuneGrid* duneGrid = Dune::GmshReader<DefaultDuneGrid>::read(fileName, verbose, insertBoundarySegments);
    return std::auto_ptr<Grid>(new DefaultGrid(duneGrid, true)); // true -> owns Dune grid
}

std::auto_ptr<Grid> importGmshGrid(const GridParameters& params,
                                   const std::string& fileName,
                                   std::vector<int>& boundaryId2PhysicalEntity,
                                   std::vector<int>& elementIndex2PhysicalEntity,
                                   bool verbose=true, bool insertBoundarySegments=true)
{
    // Check arguments
    if (params.topology != GridParameters::TRIANGULAR)
        throw std::invalid_argument("GridFactory::importGmshGrid(): unsupported grid topology");

    DefaultDuneGrid* duneGrid = Dune::GmshReader<DefaultDuneGrid>::read(
                fileName, boundaryId2PhysicalEntity, elementIndex2PhysicalEntity, verbose, insertBoundarySegments);
    return std::auto_ptr<Grid>(new DefaultGrid(duneGrid, true)); // true -> owns Dune grid
}

} // namespace Bempp
