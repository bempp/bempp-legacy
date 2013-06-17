%{
#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"
#include <algorithm>

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

#define shared_ptr boost::shared_ptr
%include "grid_factory_docstrings.i"
#undef shared_ptr

namespace Bempp
{

%ignore GridParameters;

%extend GridFactory
{


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

    static boost::shared_ptr<Bempp::Grid> createStructuredGrid(
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
    static boost::shared_ptr<Bempp::Grid> importGmshGrid(
            const std::string& topology, const std::string& fileName,
            bool verbose=false, bool insertBoundarySegments=false) {
        Bempp::GridParameters params;
        makeGridParameters(params, topology);
        return Bempp::GridFactory::importGmshGrid(params, fileName, verbose, insertBoundarySegments);
    }
    %ignore importGmshGrid;

    %apply const arma::Mat<double>& IN_MAT {
        const arma::Mat<double>& vertices
    };
    %apply const arma::Mat<int>& IN_MAT {
        const arma::Mat<int>& elementCorners
    };
    // Here we circumvent the problem that Python integers are 64-bit
    // on 64-bit system, while C++ integers may in this case be 32-bit
    static boost::shared_ptr<Bempp::Grid> createGridFromConnectivityArrays(
            const std::string& topology,
            const arma::Mat<double>& vertices,
            const arma::Mat<long>& elementCorners) {
        Bempp::GridParameters params;
        makeGridParameters(params, topology);
        arma::Mat<int> elementCornersInt(elementCorners.n_rows,
                                         elementCorners.n_cols);
        std::copy(elementCorners.begin(), elementCorners.end(),
                  elementCornersInt.begin());
        return Bempp::GridFactory::createGridFromConnectivityArrays(
            params, vertices, elementCornersInt);
    }
    %clear const arma::Mat<double>& vertices;
    %clear const arma::Mat<int>& elementCorners;
    %ignore createGridFromConnectivityArrays;
}

} // namespace Bempp
%include "grid/grid_factory.hpp"

