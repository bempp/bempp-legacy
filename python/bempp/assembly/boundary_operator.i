%{
#include "assembly/boundary_operator.hpp"
%}

// TODO
// %include "boundary_operator_docstrings.i"

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

#define shared_ptr boost::shared_ptr
%include "assembly/boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(adjoint);
}
