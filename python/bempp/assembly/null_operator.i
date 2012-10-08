%{
#include "assembly/null_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") nullOperator;
}

#define shared_ptr boost::shared_ptr
%include "assembly/null_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(nullOperator);
}
