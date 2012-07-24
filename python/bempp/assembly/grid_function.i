%{
#include "assembly/grid_function.hpp"
#include "assembly/surface_normal_dependent_function.hpp"
#include "assembly/surface_normal_independent_function.hpp"
%}

%newobject gridFunctionFromPythonSurfaceNormalIndependentFunctor;

%inline %{

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
gridFunctionFromPythonSurfaceNormalIndependentFunctor(
    const boost::shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const PythonSurfaceNormalIndependentFunctor<ResultType>& functor)
{
    return new GridFunction<BasisFunctionType, ResultType>(
        context, space, dualSpace,
        surfaceNormalIndependentFunction(functor));
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
gridFunctionFromPythonSurfaceNormalDependentFunctor(
    const boost::shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const PythonSurfaceNormalDependentFunctor<ResultType>& functor)
{
    return new GridFunction<BasisFunctionType, ResultType>(
        context, space, dualSpace,
        surfaceNormalDependentFunction(functor));
}

} // namespace Bempp

%}

namespace Bempp
{

%typemap(in) (VtkWriter::DataType dataType)
{
    if (!PyString_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "in method '$symname', argument $argnum: expected a string"); 
        SWIG_fail;
    }
    const std::string s(PyString_AsString($input));
    if (s == "cell_data")
        $1 = Bempp::VtkWriter::CELL_DATA;
    else if (s == "vertex_data")
        $1 = Bempp::VtkWriter::VERTEX_DATA;
    else
    {
        PyErr_SetString(PyExc_ValueError, "in method '$symname', argument $argnum: expected one of 'ascii', 'base64', 'appendedraw' or 'appendedbase64'");        
        SWIG_fail;
    }
}


BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);

%extend GridFunction
{
    %ignore GridFunction;

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

    %apply arma::Mat< float >& ARGOUT_MAT {
        arma::Mat< float >& result_
    };

    %apply arma::Mat< double >& ARGOUT_MAT {
        arma::Mat< double >& result_
    };

    %apply arma::Mat< std::complex<float> >& ARGOUT_MAT {
        arma::Mat<std::complex<flaot> >& result_
    };

    %apply arma::Mat< std::complex<double> >& ARGOUT_MAT {
        arma::Mat<std::complex<double> >& result_
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

    %feature("compactdefaultargs") exportToVtk;

}

%ignore gridFunctionFromFiberFunction;

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/grid_function.hpp";
#undef shared_ptr

namespace Bempp
{

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction)

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(GridFunction);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(gridFunctionFromPythonSurfaceNormalIndependentFunctor);

%clear arma::Col<float>& col_out;
%clear arma::Col<double>& col_out;
%clear arma::Col<std::complex<float> >& col_out;
%clear arma::Col<std::complex<float> >& col_out;

%clear arma::Mat<float>& result_;
%clear arma::Mat<double>& result_;
%clear arma::Mat<std::complex<float> >& result_;
%clear arma::Mat<std::complex<float> >& result_;


} // Namespace Bempp

%pythoncode %{

# functorType: "SurfaceNormalDependentFunctor" or "SurfaceNormalIndependentFunctor"
def _gridFunctionFromFunctor(
        functorType,
        context, space, dualSpace, function,
        argumentDimension, resultDimension):
    basisFunctionType = checkType(context.basisFunctionType())
    resultType = checkType(context.resultType())
    if (basisFunctionType != space.basisFunctionType() or
            basisFunctionType != dualSpace.basisFunctionType()):
        raise TypeError("BasisFunctionType of context, space and dualSpace must be the same")

    functor = constructObjectTemplatedOnValue(
        "Python" + functorType,
        resultType, function, argumentDimension, resultDimension)
    result = constructObjectTemplatedOnBasisAndResult(
        "gridFunctionFromPython" + functorType,
        basisFunctionType, resultType,
        context, space, dualSpace, functor)
    result._context = context
    result._space = space
    result._dualSpace = dualSpace
    return result

def gridFunctionFromSurfaceNormalDependentFunction(
        context, space, dualSpace, function,
        argumentDimension=3, resultDimension=1):
    return _gridFunctionFromFunctor(
        "SurfaceNormalDependentFunctor",
        context, space, dualSpace, function, 
        argumentDimension, resultDimension)

def gridFunctionFromSurfaceNormalIndependentFunction(
        context, space, dualSpace, function,
        argumentDimension=3, resultDimension=1):
    return _gridFunctionFromFunctor(
        "SurfaceNormalIndependentFunctor",
        context, space, dualSpace, function,
        argumentDimension, resultDimension)

%}
