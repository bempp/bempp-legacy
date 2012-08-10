%{
#include "linalg/default_iterative_solver.hpp"
#include "linalg/blocked_solution.hpp" // temp
#include "assembly/blocked_boundary_operator.hpp" // temp
#include "linalg/belos_solver_wrapper.hpp"
%}

// TODO
// %include "default_iterative_solver_docstrings.i"

#ifdef WITH_TRILINOS

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);

%extend DefaultIterativeSolver
{


    // add later, when we figure out how to deal with RCPs
    %ignore addPreconditioner;
}

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

%pythoncode %{

def defaultIterativeSolver(boundaryOperator):
    """Construct the default linear solver.

    This solver lets you solve the equation A f = g for the function
    f, with A being the boundary operator passed via the boundaryOperator
    argument and g a grid function supplied to the solve() method.
    """
    basisFunctionType = boundaryOperator.basisFunctionType()
    resultType = boundaryOperator.resultType()
    result = constructObjectTemplatedOnBasisAndResult(
        "DefaultIterativeSolver", basisFunctionType, resultType,
        boundaryOperator)
    result._boundaryOperator = boundaryOperator
    return result

%}

#endif
