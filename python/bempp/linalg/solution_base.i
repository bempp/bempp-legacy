%{
#include "linalg/solution_base.hpp"
%}

// TODO
// %include "solution_base_docstrings.i"

#ifdef WITH_TRILINOS

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(SolutionBase);

%extend SolutionBase
{
    %ignore thyraSolveStatus;
}

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "linalg/solution_base.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(SolutionBase);
} // namespace Bempp

#endif
