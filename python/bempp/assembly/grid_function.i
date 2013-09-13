%{
#include "assembly/grid_function.hpp"
#include "assembly/surface_normal_dependent_function.hpp"
#include "assembly/surface_normal_independent_function.hpp"
%}

%newobject gridFunctionFromPythonSurfaceNormalIndependentFunctor;
%newobject gridFunctionFromPythonSurfaceNormalDependentFunctor;

namespace Bempp {

    %apply const arma::Col<float>& IN_COL {
        const arma::Col<float>& data
        };
    %apply const arma::Col<double>& IN_COL {
        const arma::Col<double>& data
        };
    %apply const arma::Col<std::complex<float> >& IN_COL {
        const arma::Col<std::complex<float> >& data
        };
    %apply const arma::Col<std::complex<double> >& IN_COL {
        const arma::Col<std::complex<double> >& data
        };

}

%inline %{

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
uninitializedGridFunction()
{
    return new GridFunction<BasisFunctionType, ResultType>;
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
gridFunctionFromPythonSurfaceNormalIndependentFunctor(
    const boost::shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const PythonSurfaceNormalIndependentFunctor<ResultType>& functor,
    typename GridFunction<BasisFunctionType, ResultType>::ConstructionMode mode)
{
    return new GridFunction<BasisFunctionType, ResultType>(
        context, space, dualSpace,
        surfaceNormalIndependentFunction(functor),
        mode);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
gridFunctionFromPythonSurfaceNormalDependentFunctor(
    const boost::shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const PythonSurfaceNormalDependentFunctor<ResultType>& functor,
    typename GridFunction<BasisFunctionType, ResultType>::ConstructionMode mode)
{
    return new GridFunction<BasisFunctionType, ResultType>(
        context, space, dualSpace,
        surfaceNormalDependentFunction(functor),
        mode);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
gridFunctionFromCoefficients(
    const boost::shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const arma::Col<ResultType>& data)
{
    return new GridFunction<BasisFunctionType, ResultType>(
        context, space, dualSpace, data,
        GridFunction<BasisFunctionType, ResultType>::COEFFICIENTS);
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>*
gridFunctionFromProjections(
    const boost::shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const boost::shared_ptr<const Space<BasisFunctionType> >& space,
    const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
    const arma::Col<ResultType>& data)
{
    return new GridFunction<BasisFunctionType, ResultType>(
        context, space, dualSpace, data, 
        GridFunction<BasisFunctionType, ResultType>::PROJECTIONS);
}


} // namespace Bempp

%}

namespace Bempp
{

%typemap(typecheck) VtkWriter::DataType
{
    $1 = PyString_Check($input);
}
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
        PyErr_SetString(PyExc_ValueError,
                        "in method '$symname', argument $argnum: "
                        "expected either 'cell_data' or 'vertex_data'");
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

    %apply arma::Mat< float >& ARGOUT_MAT {
        arma::Mat< float >& points
    };

    %apply arma::Mat< double >& ARGOUT_MAT {
        arma::Mat< double >& points
    };
  
    %apply arma::Mat< float >& ARGOUT_MAT {
        arma::Mat< float >& values
    };

    %apply arma::Mat< double >& ARGOUT_MAT {
        arma::Mat< double >& values
    };

    %apply arma::Mat< std::complex<float> >& ARGOUT_MAT {
        arma::Mat<std::complex<flaot> >& values
    };

    %apply arma::Mat< std::complex<double> >& ARGOUT_MAT {
        arma::Mat<std::complex<double> >& values
    };

    %ignore evaluateAtSpecialPoints;
    void _evaluateAtSpecialPoints(
            VtkWriter::DataType dataType,
            arma::Mat<CoordinateType>& points, arma::Mat<ResultType>& values) const {
        $self->evaluateAtSpecialPoints(dataType, points, values); 
    }

    %pythoncode %{
    def evaluateAtSpecialPoints(self, dataType, returnPoints=False):
       points, values = self._evaluateAtSpecialPoints(dataType)
       if returnPoints:
          return points, values
       else:
          return values
    %}
    
    void coefficients(arma::Col<ResultType>& col_out)
    {
        col_out = $self->coefficients();
    }

    void projections(arma::Col<ResultType>& col_out)
    {
        col_out = $self->projections();
    }

    void projections(const boost::shared_ptr<const Space<BasisFunctionType> >& dualSpace,
                     arma::Col<ResultType>& col_out)
    {
        col_out = $self->projections(dualSpace);
    }

    %ignore coefficients;
    %ignore projections;

    GridFunction<BasisFunctionType, ResultType> __pos__()
    {
        return +(*$self);
    }

    GridFunction<BasisFunctionType, ResultType> __neg__()
    {
        return -(*$self);
    }

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

    %pythoncode %{
      def plot(self):
          "Visualize the GridFunction."

          from visualization import plotGridFunction
          plotGridFunction(self)
	    %}
   
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
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(uninitializedGridFunction);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(gridFunctionFromPythonSurfaceNormalIndependentFunctor);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(gridFunctionFromPythonSurfaceNormalDependentFunctor);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(gridFunctionFromCoefficients);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(gridFunctionFromProjections);
%feature("compactdefaultargs") exportToVtk;
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(exportToVtk);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(exportToGmsh);



%clear arma::Col<float>& col_out;
%clear arma::Col<double>& col_out;
%clear arma::Col<std::complex<float> >& col_out;
%clear arma::Col<std::complex<double> >& col_out;

%clear arma::Mat<float>& result_;
%clear arma::Mat<double>& result_;
%clear arma::Mat<std::complex<float> >& result_;
%clear arma::Mat<std::complex<double> >& result_;

%clear arma::Mat<float>& points;
%clear arma::Mat<double>& points;
%clear arma::Mat<float>& values;
%clear arma::Mat<double>& values;
%clear arma::Mat<std::complex<float> >& values;
%clear arma::Mat<std::complex<double> >& values;
 
%clear arma::Mat<float>& data;
%clear arma::Mat<double>& data;
%clear arma::Mat<std::complex<float> >& data;
%clear arma::Mat<std::complex<double> >& data;


} // Namespace Bempp
