%{
#include "assembly/grid_function.hpp"
%}

%define BEMPP_DECLARE_GRID_FUNCTION(CLASS1,CLASS2,NPY1,NPY2)

namespace Bempp
{

%extend GridFunction< CLASS1 , CLASS2 >
 {

    %ignore coefficients;
    %ignore projections;
    %ignore setCoefficients;
    %ignore setProjections;
    %ignore grid;
    %ignore codomainDimension;
    %ignore space;

    %ignore basis;
    %ignore getLocalCoefficients;
    %ignore exportToVtk;

 }


 %ignore gridFunctionFromFiberFunction;
 %ignore gridFunctionFromSurfaceNormalDependentFunctor;

} // Bempp

%apply const arma::Col< CLASS2 >& IN_COL { const std::arma::Col< CLASS2 >& coefficients };
%apply const arma::Col< CLASS2 >& IN_COL { const std::arma::Col< CLASS2 >& projections };


%include "assembly/grid_function.hpp";

namespace Bempp
{
  %template(GridFunction_## NPY1 ##_## NPY2 ) GridFunction< CLASS1 , CLASS2 >;
  %template(gridFunctionFromSurfaceNormalIndependentFunctor_## NPY1 ##_## NPY2) 
     gridFunctionFromSurfaceNormalIndependentFunctor< CLASS1 , CLASS2 >;
}

%clear const std::arma::Col< CLASS2 >& coefficients;
%clear const std::arma::Col< CLASS2 >& projections;
%enddef // BEMPP_DECLARE_GRID_FUNCTION

ITERATE_OVER_BASIS_RESULT_TYPES(BEMPP_DECLARE_GRID_FUNCTION)

/*

%define GRID_FUNCTION_APPLY_MAPS_ENABLE(TYPE)
%enddef

%define GRID_FUNCTION_APPLY_MAPS_DISABLE(TYPE)
%enddef

ITERATE_OVER_TYPES(GRID_FUNCTION_APPLY_MAPS_ENABLE)

%include "assembly/grid_function.hpp"


ITERATE_OVER_TYPES(GRID_FUNCTION_APPLY_MAPS_DISABLE)

*/
