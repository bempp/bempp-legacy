%{
#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"
#include "./grid/grid_parameters_converter.hpp"
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
    
    %feature("compactdefaultargs") importGmshGrid;
    static std::auto_ptr<Bempp::Grid> importGmshGrid(
            const std::string& topology, const std::string& fileName, 
            bool verbose=true, bool insertBoundarySegments=true) {
        Bempp::GridParameters params;
        makeGridParameters(params, topology);
        return Bempp::GridFactory::importGmshGrid(params, fileName, verbose, insertBoundarySegments);
    }
    %ignore importGmshGrid;
}

} // namespace Bempp
%include "grid/grid_factory.hpp"

