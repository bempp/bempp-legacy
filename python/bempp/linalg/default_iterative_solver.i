%{
#include "linalg/default_iterative_solver.hpp"
#include "linalg/blocked_solution.hpp" // temp
#include "linalg/belos_solver_wrapper.hpp"
%}

// TODO
// %include "default_iterative_solver_docstrings.i"

#ifdef WITH_TRILINOS

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);

%define BEMPP_EXTEND_DEFAULT_ITERATIVE_SOLVER(BASIS, RESULT, PYBASIS, PYRESULT)
    %extend DefaultIterativeSolver< BASIS, RESULT >
    {
        %pythonprepend setPreconditioner %{
            print ("DefaultIterativeSolver.setPreconditioner() is deprecated. "
                   "Pass the preconditioner directly to the initializeSolver() "
                   "method.")
        %}
    }
%enddef
BEMPP_ITERATE_OVER_BASIS_AND_RESULT_TYPES(BEMPP_EXTEND_DEFAULT_ITERATIVE_SOLVER);

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(
    double tol, int maxIterationCount = 1000);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(
    double tol, int maxIterationCount = 1000);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "linalg/default_iterative_solver.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);
} // namespace Bempp

#endif // WITH_TRILINOS
