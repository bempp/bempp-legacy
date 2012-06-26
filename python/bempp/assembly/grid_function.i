%{
#include "assembly/grid_function.hpp"
%}

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

//%define BEMPP_DECLARE_GRID_FUNCTION(CLASS1,CLASS2,NPY1,NPY2)

%extend GridFunction
{
    %ignore setCoefficients;
    %ignore setProjections;
    %ignore codomainDimension;
    %ignore space;

    %ignore basis;
    %ignore getLocalCoefficients;

    %apply arma::Col<float>& ARGOUT_COL {
        arma::Col<float>& col_out
    };
    %apply arma::Col<double>& ARGOUT_COL {
        arma::Col<double>& col_out
    };
    %apply arma::Col<std::complex<float> >& ARGOUT_COL {
        arma::Col<std::complex<float> >& col_out
    };
    %apply arma::Col<std::complex<double> >& ARGOUT_COL {
        arma::Col<std::complex<double> >& col_out
    };

    void coefficients(arma::Col<ResultType>& col_out)
    {
        col_out = $self->coefficients();
    }

    void projections(arma::Col<ResultType>& col_out)
    {
        col_out = $self->projections();
    }

    %ignore coefficients;
    %ignore projections;
}

%ignore gridFunctionFromFiberFunction;

} // namespace Bempp

%include "assembly/grid_function.hpp";

namespace Bempp
{

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    gridFunctionFromSurfaceNormalIndependentFunctor);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    gridFunctionFromSurfaceNormalDependentFunctor);

%clear arma::Col<float>& col_out;
%clear arma::Col<double>& col_out;
%clear arma::Col<std::complex<float> >& col_out;
%clear arma::Col<std::complex<float> >& col_out;

} // namespace Bempp

%pythoncode %{

def gridFunctionFromSurfaceNormalDependentFunctor(space, functor, factory,
        assemblyOptions):
    basisFunctionType = checkType(space.basisFunctionType())
    resultType = checkType(functor.valueType())
    result = constructObjectTemplatedOnBasisAndResult(
        "gridFunctionFromSurfaceNormalDependentFunctor",
        basisFunctionType, resultType,
        space, functor, factory, assemblyOptions)
    result._space = space
    return result

def gridFunctionFromSurfaceNormalIndependentFunctor(space, functor, factory,
        assemblyOptions, basisFunctionType='float64'):
    basisFunctionType = checkType(space.basisFunctionType())
    resultType = checkType(functor.valueType())
    result = constructObjectTemplatedOnBasisAndResult(
        "gridFunctionFromSurfaceNormalIndependentFunctor",
        basisFunctionType, resultType,
        space, functor, factory, assemblyOptions)
    result._space = space
    return result

%}
