%{
#include "linalg/solver.hpp"
%}

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
}

#define shared_ptr boost::shared_ptr
%include "linalg/solver.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
}
