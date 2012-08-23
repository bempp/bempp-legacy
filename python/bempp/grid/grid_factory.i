%{
#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"

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
    else if (topology == "hybrid_2d")
        params.topology = Bempp::GridParameters::HYBRID_2D;
    else if (topology == "tetrahedral")
        params.topology = Bempp::GridParameters::TETRAHEDRAL;
    else
        throw std::runtime_error("Invalid grid topology requested");
}

} // namespace Bempp
%}

%include "grid_factory_docstrings.i"

namespace Bempp 
{
  
%ignore GridParameters;
  
%extend GridFactory 
{

    // Note: no %newobject directives are necessary; auto_ptr's typemaps ensure
    // that ownership is passed to Python.
    %apply const arma::Col<ctype>& IN_COL {
        const arma::Col<ctype>& lowerLeft, 
        const arma::Col<ctype>& upperRight 
    };
    %apply const arma::Col<unsigned int>& IN_COL {
        const arma::Col<unsigned int>& nElements
    };

  //    %pythonappend createStructuredGrid %{
  //    val.topology = args[0]
  //	%}

    static std::auto_ptr<Bempp::Grid> createStructuredGrid(
            const std::string& topology,
            const arma::Col<Bempp::ctype>& lowerLeft, 
            const arma::Col<Bempp::ctype>& upperRight,
            const arma::Col<unsigned int>& nElements) {
        Bempp::GridParameters params;
        makeGridParameters(params, topology);
        return Bempp::GridFactory::createStructuredGrid(params, lowerLeft, upperRight, nElements);
    }
    %clear const arma::Col<ctype>& lowerLeft; 
    %clear const arma::Col<ctype>& upperRight;
    %clear const arma::Col<unsigned int>& nElements;
    %ignore createStructuredGrid;

    //%pythonappend importGmshGrid %{
    //  val.topology = args[0]
    //	%}
    
    %feature("compactdefaultargs") importGmshGrid;
    static std::auto_ptr<Bempp::Grid> importGmshGrid(
            const std::string& topology, const std::string& fileName, 
            bool verbose=true, bool insertBoundarySegments=false) {
        Bempp::GridParameters params;
        makeGridParameters(params, topology);
        return Bempp::GridFactory::importGmshGrid(params, fileName, verbose, insertBoundarySegments);
    }
    %ignore importGmshGrid;
}

} // namespace Bempp
%include "grid/grid_factory.hpp"

