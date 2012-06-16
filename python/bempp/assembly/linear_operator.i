%{
#include "assembly/linear_operator.hpp"
#include "assembly/linear_operator_superposition.hpp"
#include <complex>
%}

// TODO
// %include "linear_operator_docstrings.i"

namespace Bempp
{

BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);

template <typename BasisFunctionType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend LinearOperator
{
    // this function is only for internal use
    %ignore constituentOperators;

    // this function is only for internal use
    %ignore constituentOperatorWeights;

    LinearOperatorSuperposition<BasisFunctionType, ResultType> __add__(
        const LinearOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self + other;
    }

    LinearOperatorSuperposition<BasisFunctionType, ResultType> __sub__(
        const LinearOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self - other;
    }

    LinearOperatorSuperposition<BasisFunctionType, ResultType> __mul__(
        ResultType other)
    {
        return *$self * other;
    }

    LinearOperatorSuperposition<BasisFunctionType, ResultType> __rmul__(
        ResultType other)
    {
        return *$self * other;
    }

    LinearOperatorSuperposition<BasisFunctionType, ResultType> __div__(
        ResultType other)
    {
        return *$self / other;
    }
}

BEMPP_PYTHON_EXTEND_INTERFACE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);

} // namespace Bempp

%include "assembly/linear_operator.hpp"

namespace Bempp
{
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);
}
