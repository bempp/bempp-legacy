%{
#include "assembly/elementary_integral_operator.hpp"
%}

// TODO
// %include "elementary_integral_operator_docstrings.i"

%include "assembly/elementary_integral_operator.hpp"

namespace Bempp
{
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ElementaryIntegralOperator);
}
