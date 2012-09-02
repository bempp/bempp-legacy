%{
#include "assembly/abstract_boundary_operator.hpp"
#include <complex>
%}

// TODO
// %include "abstract_boundary_operator_docstrings.i"

BEMPP_DECLARE_SHARED_PTR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    Bempp::AbstractBoundaryOperator);

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename BasisFunctionType> class Space;
class AssemblyOptions;
class Symmetry;

%extend AbstractBoundaryOperator
{
    %ignore id;
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/abstract_boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperator);
}
