%{
#include "assembly/identity_operator.hpp"
%}

#define shared_ptr boost::shared_ptr
%include "assembly/identity_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(identityOperator);
}

%pythoncode %{

def identityOperator(context, domain, range, dualToRange):
    """Construct an identity operator."""
    return _constructOperator(
        "identityOperator", context, domain, range, dualToRange)

%}
