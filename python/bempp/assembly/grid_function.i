%{
#include "assembly/grid_function.hpp"
%}

%inline %{

namespace Bempp
{
template <typename BasisFunctionType, typename Functor>
GridFunction<BasisFunctionType, typename Functor::ValueType>
inline gridFunctionFromSurfaceNormalIndependentFunctorX(
    const boost::shared_ptr<const Context<BasisFunctionType, typename Functor::ValueType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const Functor& functor)
{
    gridFunctionFromSurfaceNormalIndependentFunctor(context, space, dualSpace, functor);
}

template <typename BasisFunctionType, typename Functor>
GridFunction<BasisFunctionType, typename Functor::ValueType>
inline gridFunctionFromSurfaceNormalDependentFunctorX(
    const boost::shared_ptr<const Context<BasisFunctionType, typename Functor::ValueType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const Functor& functor)
{
    gridFunctionFromSurfaceNormalDependentFunctor(context, space, dualSpace, functor);
}

} // namespace Bempp

%}

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

%define BEMPP_INSTANTIATE_GRID_FUNCTION_FROM_FUNCTOR(FUNCTOR)
    %template(pythonGridFunctionFrom ## FUNCTOR ## _float32_float32)
        gridFunctionFrom ## FUNCTOR ## X<float, Python ## FUNCTOR<float> >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _float32_complex64)
        gridFunctionFrom ## FUNCTOR ## X<float, Python ## FUNCTOR<std::complex<float> > >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _complex64_complex64)
        gridFunctionFrom ## FUNCTOR ## X<std::complex<float>, Python ## FUNCTOR<std::complex<float> > >;

    %template(pythonGridFunctionFrom ## FUNCTOR ## _float64_float64)
        gridFunctionFrom ## FUNCTOR ## X<double, Python ## FUNCTOR<double> >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _float64_complex128)
        gridFunctionFrom ## FUNCTOR ## X<double, Python ## FUNCTOR<std::complex<double> > >;
    %template(pythonGridFunctionFrom ## FUNCTOR ## _complex128_complex128)
        gridFunctionFrom ## FUNCTOR ## X<std::complex<double>, Python ## FUNCTOR<std::complex<double> > >
%enddef

%extend GridFunction
{
    %ignore GridFunction;

    boost::shared_ptr<const Space<BasisFunctionType> > space() const
    {
        return $self->space();
    }
    %ignore space;

    boost::shared_ptr<const Space<BasisFunctionType> > dualSpace() const
    {
        return $self->dualSpace();
    }
    %ignore dualSpace;

    boost::shared_ptr<const Context<BasisFunctionType, ResultType> > context() const
    {
        return $self->context();
    }
    %ignore context;

    %ignore setCoefficients;
    %ignore setProjections;
    %ignore codomainDimension;

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

    GridFunction<BasisFunctionType, ResultType> __add__(
        const GridFunction<BasisFunctionType, ResultType>& other)
    {
        return *$self + other;
    }

    GridFunction<BasisFunctionType, ResultType> __sub__(
        const GridFunction<BasisFunctionType, ResultType>& other)
    {
        return *$self - other;
    }

    GridFunction<BasisFunctionType, ResultType> __mul__(
        ResultType other)
    {
        return *$self * other;
    }

    GridFunction<BasisFunctionType, ResultType> __rmul__(
        ResultType other)
    {
        return *$self * other;
    }

    GridFunction<BasisFunctionType, ResultType> __div__(
        ResultType other)
    {
        return *$self / other;
    }
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

%clear arma::Col<float>& col_out;
%clear arma::Col<double>& col_out;
%clear arma::Col<std::complex<float> >& col_out;
%clear arma::Col<std::complex<float> >& col_out;

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
