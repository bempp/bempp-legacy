%{
#include "assembly/identity_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") identityOperator;
}

#define shared_ptr boost::shared_ptr
%include "assembly/identity_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(identityOperator);
}
