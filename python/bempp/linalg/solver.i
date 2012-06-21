%{
#include "linalg/solver.hpp"
%}

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);
}

%include "linalg/solver.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Solver);
}
