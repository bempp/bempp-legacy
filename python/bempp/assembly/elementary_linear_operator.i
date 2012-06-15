%{
#include "assembly/elementary_linear_operator.hpp"
%}

// TODO
// %include "elementary_linear_operator_docstrings.i"

namespace Bempp
{

BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ElementaryLinearOperator);

%extend ElementaryLinearOperator
{
    // this function is only for internal use
    %ignore makeAssembler;

    // this function is only for internal use
    %ignore assembleDetachedWeakFormInternal;
}

} // namespace Bempp

%include "assembly/elementary_linear_operator.hpp"

namespace Bempp
{
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ElementaryLinearOperator);
}
