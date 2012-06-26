%{
#include "assembly/scaled_linear_operator.hpp"
%}

// TODO
// %include "scaled_linear_operator_docstrings.i"

%include "assembly/scaled_linear_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(ScaledLinearOperator);
}
