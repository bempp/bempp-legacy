%{
#include "assembly/blocked_operator_structure.hpp"
  %}

namespace Bempp
{
  BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedOperatorStructure);
  BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedOperatorStructure);
}

%include "assembly/blocked_operator_structure.hpp"

namespace Bempp
{
  BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(BlockedOperatorStructure);
}

