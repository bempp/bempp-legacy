%{
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/grid_parameters.hpp"
#include "common/shared_ptr.hpp"
%}

%include "grid_docstrings.i"

%shared_ptr(Bempp::Grid);

%include "std_vector.i"
namespace std
{
%template(vector_bool) vector<bool>;
}

namespace Bempp
{

%extend Grid
{
    %pythonappend leafView %{
        val._parentGrid = self
    %}

    %pythonappend levelView %{
        val._parentGrid = self
    %}

    %pythonappend globalIdSet %{
        val._parentGrid = self
    %}

    %pythoncode %{
      def plot(self):
          "Visualize the Grid."

          from visualization import plotGrid
          plotGrid(self)
        %}


    // this function is only for internal use
    %ignore elementGeometryFactory;
}

%apply const arma::Mat<float>& IN_MAT {
    const arma::Mat<float >& points
};
%apply const arma::Mat<double>& IN_MAT {
    const arma::Mat<double >& points
};

%ignore areInside(const Grid& grid, const arma::Mat<float>& points);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "grid/grid.hpp"
#undef shared_ptr


%clear arma::Mat<float>& points;
%clear arma::Mat<double>& points;
