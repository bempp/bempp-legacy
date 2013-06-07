%{
#include "linalg/default_direct_solver.hpp"
#include "linalg/blocked_solution.hpp" // temp
%}

// TODO
// %include "default_direct_solver_docstrings.i"

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultDirectSolver);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "linalg/default_direct_solver.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(DefaultDirectSolver);
} // namespace Bempp
