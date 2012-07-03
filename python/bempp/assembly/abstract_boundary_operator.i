%{
#include "assembly/abstract_boundary_operator.hpp"
#include "assembly/abstract_boundary_operator_sum.hpp"
#include <complex>
%}

// TODO
// %include "abstract_boundary_operator_docstrings.i"

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend AbstractBoundaryOperator
{
    %ignore clone;

    // SWIG doesn't handle properly the using boost::shared_ptr declaration
    // in BEM++'s code, therefore this wrapper returning a fully qualified
    // shared_ptr is needed
    boost::shared_ptr<const DiscreteBoundaryOperator<ResultType> > weakForm() const
    {
        return $self->weakForm();
    }
    %ignore weakForm;

    AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> __add__(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self + other;
    }

    AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> __sub__(
        const AbstractBoundaryOperator<BasisFunctionType, ResultType>& other)
    {
        return *$self - other;
    }

    ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> __mul__(
        ResultType other)
    {
        return *$self * other;
    }

    GridFunction<BasisFunctionType, ResultType> __mul__(
        const GridFunction<BasisFunctionType, ResultType>& other)
    {
        return *$self * other;
    }

    ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> __rmul__(
        ResultType other)
    {
        return *$self * other;
    }

    ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> __div__(
        ResultType other)
    {
        return *$self / other;
    }
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

} // namespace Bempp

%include "assembly/abstract_boundary_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);
}
