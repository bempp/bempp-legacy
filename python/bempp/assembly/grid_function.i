%{
#include "assembly/grid_function.hpp"
%}

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

%define BEMPP_INSTANTIATE_GRID_FUNCTION_FROM_FUNCTOR(FUNCTOR)
    %template(pythonGridFunctionFrom ## FUNCTOR ## _float32_float32)
        gridFunctionFrom ## FUNCTOR<float, Python ## FUNCTOR<float> >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _float32_complex64)
        gridFunctionFrom ## FUNCTOR<float, Python ## FUNCTOR<std::complex<float> > >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _complex64_complex64)
        gridFunctionFrom ## FUNCTOR<std::complex<float>, Python ## FUNCTOR<std::complex<float> > >;

    %template(pythonGridFunctionFrom ## FUNCTOR ## _float64_float64)
        gridFunctionFrom ## FUNCTOR<double, Python ## FUNCTOR<double> >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _float64_complex128)
        gridFunctionFrom ## FUNCTOR<double, Python ## FUNCTOR<std::complex<double> > >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _complex128_complex128)
        gridFunctionFrom ## FUNCTOR<std::complex<double>, Python ## FUNCTOR<std::complex<double> > >
%enddef

%extend GridFunction
{
    %ignore coefficients;
    %ignore projections;
    %ignore setCoefficients;
    %ignore setProjections;
    %ignore codomainDimension;

    %ignore basis;
    %ignore getLocalCoefficients;
}

%ignore gridFunctionFromFiberFunction;

} // namespace Bempp

%include "assembly/grid_function.hpp";

namespace Bempp
{

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction)

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);
BEMPP_INSTANTIATE_GRID_FUNCTION_FROM_FUNCTOR(
    SurfaceNormalIndependentFunctor);
BEMPP_INSTANTIATE_GRID_FUNCTION_FROM_FUNCTOR(
    SurfaceNormalDependentFunctor);

} // namespace Bempp

%pythoncode %{

# functorType: "SurfaceNormalDependentFunctor" or "SurfaceNormalIndependentFunctor"
def _gridFunctionFromFunctor(
        functorType,
        space, dualSpace, function, factory, assemblyOptions,
        valueType, argumentDimension, resultDimension):
    basisFunctionType = checkType(space.basisFunctionType())
    resultType = checkType(valueType)
    functor = constructObjectTemplatedOnValue(
        "Python" + functorType,
        valueType, function, argumentDimension, resultDimension)
    result = constructObjectTemplatedOnBasisAndResult(
        "pythonGridFunctionFrom" + functorType,
        basisFunctionType, resultType,
        space, dualSpace, functor, factory, assemblyOptions)
    result._space = space
    result._dualSpace = dualSpace
    return result

def gridFunctionFromSurfaceNormalDependentFunction(
        space, dualSpace, function, factory, assemblyOptions,
        valueType='float64', argumentDimension=3, resultDimension=1):
    return _gridFunctionFromFunctor(
        "SurfaceNormalDependentFunctor",
        space, dualSpace, function, factory, assemblyOptions,
        valueType, argumentDimension, resultDimension)

def gridFunctionFromSurfaceNormalIndependentFunction(
        space, dualSpace, function, factory, assemblyOptions,
        valueType='float64', argumentDimension=3, resultDimension=1):
    return _gridFunctionFromFunctor(
        "SurfaceNormalIndependentFunctor",
        space, dualSpace, function, factory, assemblyOptions,
        valueType, argumentDimension, resultDimension)

%}
