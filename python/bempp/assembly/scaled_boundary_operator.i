%{
#include "assembly/scaled_boundary_operator.hpp"
%}

// TODO
// %include "scaled_boundary_operator_docstrings.i"

%include "assembly/scaled_boundary_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(ScaledBoundaryOperator);
}
