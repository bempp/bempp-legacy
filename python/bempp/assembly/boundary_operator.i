%{
#include "assembly/boundary_operator.hpp"
#include <complex>
%}

// TODO
// %include "boundary_operator_docstrings.i"

%shared_ptr(Bempp::BoundaryOperator<float,float>);
%shared_ptr(Bempp::BoundaryOperator<float,std::complex<float> >);
%shared_ptr(Bempp::BoundaryOperator<double,double >);
%shared_ptr(Bempp::BoundaryOperator<double,std::complex<double> >);
%shared_ptr(Bempp::BoundaryOperator<std::complex<float>,std::complex<float> >);
%shared_ptr(Bempp::BoundaryOperator<std::complex<double>,std::complex<double> >);

#define shared_ptr boost::shared_ptr
namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend BoundaryOperator
{
    %ignore BoundaryOperator;

    BoundaryOperator<BasisFunctionType, ResultType> __pos__()
    {
        return +(*$self);
    }

    BoundaryOperator<BasisFunctionType, ResultType> __neg__()
    {
        return -(*$self);
    }

    BoundaryOperator<BasisFunctionType, ResultType> __add__(
        const BoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self + other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __sub__(
        const BoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self - other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __mul__(
        ResultType other)
    {
        return *$self * other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __mul__(
        const BoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self * other;
    }


    GridFunction<BasisFunctionType, ResultType> __mul__(
        const GridFunction<BasisFunctionType, ResultType>& other)
    {
        return *$self * other;
    }


    BoundaryOperator<BasisFunctionType, ResultType> __rmul__(
        ResultType other)
    {
        return *$self * other;
    }

    BoundaryOperator<BasisFunctionType, ResultType> __div__(
        ResultType other)
    {
        return *$self / other;
    }
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);

} // namespace Bempp

%include "assembly/boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(adjoint);
}
