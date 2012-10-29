%{
#include "assembly/identity_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") identityOperator;

// Do not emit warnings about ignored assignment operators.
%warnfilter(362) IdentityOperator::operator=;
}

#define shared_ptr boost::shared_ptr
%include "assembly/identity_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(identityOperator);
}
