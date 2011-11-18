#ifndef bempp_lib_grid_grid_factory_hpp
#define bempp_lib_grid_grid_factory_hpp

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
