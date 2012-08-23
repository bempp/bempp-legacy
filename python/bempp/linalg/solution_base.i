%{
#include "linalg/solution_base.hpp"
%}

// TODO
// %include "solution_base_docstrings.i"

#ifdef WITH_TRILINOS

namespace Bempp
{

%typemap(out) SolutionStatus::Status
{
    if ($1 == Bempp::SolutionStatus::CONVERGED)
        $result = PyString_FromString("converged");
    else if ($1 == Bempp::SolutionStatus::UNCONVERGED)
        $result = PyString_FromString("unconverged");
    else if ($1 ==Bempp::SolutionStatus::UNKNOWN)
        $result = PyString_FromString("unknown");
    else
    {
        PyErr_SetString(PyExc_ValueError, "in method '$symname', unknown solver status.");
        SWIG_fail;
    }
}

}

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(SolutionBase);

%extend SolutionBase
{
    %ignore thyraSolveStatus;
    %ignore unknownTolerance;
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
