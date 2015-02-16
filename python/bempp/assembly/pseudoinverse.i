%{
#include "assembly/abstract_boundary_operator_pseudoinverse.hpp"
%}

#define shared_ptr boost::shared_ptr
%include "assembly/abstract_boundary_operator_pseudoinverse.hpp"
#undef

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(pseudoinverse);
}

