%{
#include "linalg/blocked_solution.hpp"
%}

// TODO
// %include "blocked_solution_docstrings.i"

#ifdef WITH_TRILINOS

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedSolution);

%extend BlockedSolution
{
%pythonappend gridFunction %{
    val._parentSolution = self
%}
}

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "linalg/blocked_solution.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(BlockedSolution);
} // namespace Bempp

#endif
