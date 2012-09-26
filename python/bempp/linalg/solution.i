%{
#include "linalg/solution.hpp"
%}

// TODO
// %include "solution_docstrings.i"

#ifdef WITH_TRILINOS

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solution);

%extend Solution
{

%pythonappend gridFunction %{
    val._parentSolution = self
%}

}

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "linalg/solution.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Solution);
} // namespace Bempp

#endif
