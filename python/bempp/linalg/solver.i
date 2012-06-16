%{
#include "linalg/solver.hpp"
%}

%include "linalg/solver.hpp"

namespace Bempp
{

BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);

} // namespace Bempp
