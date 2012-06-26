%{
#include "assembly/linear_operator.hpp"
#include "assembly/linear_operator_sum.hpp"
#include <complex>
%}

// TODO
// %include "linear_operator_docstrings.i"

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend LinearOperator
{
    %ignore clone;

    // SWIG doesn't handle properly the using boost::shared_ptr declaration
    // in BEM++'s code, therefore this wrapper returning a fully qualified
    // shared_ptr is needed
    boost::shared_ptr<const DiscreteLinearOperator<ResultType> > weakForm() const
    {
        return $self->weakForm();
    }
    %ignore weakForm;

    LinearOperatorSum<BasisFunctionType, ResultType> __add__(
        const LinearOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self + other;
    }

    LinearOperatorSum<BasisFunctionType, ResultType> __sub__(
        const LinearOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self - other;
    }

    ScaledLinearOperator<BasisFunctionType, ResultType> __mul__(
        ResultType other)
    {
        return *$self * other;
    }

    GridFunction<BasisFunctionType, ResultType> __mul__(
        const GridFunction<BasisFunctionType, ResultType>& other)
    {
        return *$self * other;
    }

    ScaledLinearOperator<BasisFunctionType, ResultType> __rmul__(
        ResultType other)
    {
        return *$self * other;
    }

    ScaledLinearOperator<BasisFunctionType, ResultType> __div__(
        ResultType other)
    {
        return *$self / other;
    }
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);

} // namespace Bempp

%include "assembly/linear_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperator);
}
