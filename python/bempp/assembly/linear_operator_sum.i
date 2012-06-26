%{
#include "assembly/linear_operator_sum.hpp"
%}

// TODO
// %include "linear_operator_sum_docstrings.i"

%include "assembly/linear_operator_sum.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperatorSum);
}
