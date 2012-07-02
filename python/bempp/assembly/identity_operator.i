%{
#include "assembly/identity_operator.hpp"
%}

// TODO
// %include "identity_operator_docstrings.i"

%include "assembly/identity_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(IdentityOperator);
}

%pythoncode %{

def identityOperator(domain, range, dualToRange, resultType=None):
    """Construct an identity operator."""
    return _constructOperator(
    "IdentityOperator", domain, range, dualToRange, resultType)

%}
