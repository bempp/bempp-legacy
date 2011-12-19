#ifndef bempp_grid_parameters_converter_hpp
#define bempp_grid_parameters_converter_hpp

#include "grid/grid_factory.hpp"

namespace Bempp 
{

inline void makeGridParameters(GridParameters& params, const std::string& topology)
{
    if (topology == "linear")
        params.topology = Bempp::GridParameters::LINEAR;
    else if (topology == "triangular")
        params.topology = Bempp::GridParameters::TRIANGULAR;
    else if (topology == "quadrilateral")
        params.topology = Bempp::GridParameters::QUADRILATERAL;
    else if (topology == "hybrid")
        params.topology = Bempp::GridParameters::HYBRID;
    else
        throw std::runtime_error("Invalid grid topology requested");
}

} // namespace Bempp

#endif
