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



 }


 %ignore gridFunctionFromFiberFunction;

} // Bempp

%include "assembly/grid_function.hpp";

namespace Bempp
{
  %template(GridFunction_## NPY1 ##_## NPY2 ) GridFunction< CLASS1 , CLASS2 >;
  %template(gridFunctionFromSurfaceNormalIndependentFunctor_## NPY1 ##_## NPY2) 
     gridFunctionFromSurfaceNormalIndependentFunctor< CLASS1 , CLASS2 >;
  %template(gridFunctionFromSurfaceNormalDependentFunctor_## NPY1 ##_## NPY2) 
     gridFunctionFromSurfaceNormalDependentFunctor< CLASS1 , CLASS2 >;
}

%enddef // BEMPP_DECLARE_GRID_FUNCTION

ITERATE_OVER_BASIS_RESULT_TYPES(BEMPP_DECLARE_GRID_FUNCTION)

%pythoncode %{

def gridFunctionFromSurfaceNormalDependentFunctor(space,functor,factory,assemblyOptions,basisType='float64',valueType='float64'):
        dtype1=checkType(basisType)
	dtype2=checkType(valueType)
        if not dtype2==functor._valueType:
                raise ValueError("ValueType of functor is "+functor._valueType+"!")
	fullName="gridFunctionFromSurfaceNormalDependentFunctor"+"_"+dtype1+"_"+dtype2
	keys=globals()
	if not fullName in keys: raise NameError("Could not find "+fullName+"!")
	tmp=keys[fullName](space,functor,factory,assemblyOptions)
        tmp._basisType=checkType(basisType)
        tmp._valueType=checkType(valueType)
        return tmp

def gridFunctionFromSurfaceNormalIndependentFunctor(space,functor,factory,assemblyOptions,basisType='float64',valueType='float64'):
        dtype1=checkType(basisType)
	dtype2=checkType(valueType)
        if not dtype2==functor._valueType:
                raise ValueError("ValueType of functor is "+functor._valueType+"!")
	fullName="gridFunctionFromSurfaceNormalIndependentFunctor"+"_"+dtype1+"_"+dtype2
	keys=globals()
	if not fullName in keys: raise NameError("Could not find "+fullName+"!")
	tmp=keys[fullName](space,functor,factory,assemblyOptions)
        tmp._basisType=dtype1
        tmp._valueType=dtype2
        return tmp

  %}
