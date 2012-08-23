%{
#include "assembly/blocked_boundary_operator.hpp"
  %}

namespace Bempp
{
  BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedBoundaryOperator);

  %extend BlockedBoundaryOperator {
    %ignore apply;
  }

  BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedBoundaryOperator);
}

#define shared_ptr boost::shared_ptr
%include "assembly/blocked_boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{
  BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(BlockedBoundaryOperator);
}
