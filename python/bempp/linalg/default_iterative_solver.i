%{
#include "linalg/default_iterative_solver.hpp"
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

    %ignore getThyraSolveStatus;
}

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol);

} // namespace Bempp

%include "linalg/default_iterative_solver.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);
} // namespace Bempp

%pythoncode %{

def defaultIterativeSolver(linearOperator, gridFunction):
    """Construct the default linear solver.

    This solver lets you solve the equation A f = g for the function
    f, with A being the linear operator passed via the linearOperator
    argument and g the function represented by the gridFunction
    argument
    """
    basisFunctionType = linearOperator.basisFunctionType()
    if (basisFunctionType != gridFunction.basisFunctionType()):
        raise TypeError("BasisFunctionType of linearOperator and "
                        "gridFunction must be the same")
    resultType = linearOperator.resultType()
    if (resultType != gridFunction.resultType()):
        raise TypeError("ResultType of linearOperator and "
                        "gridFunction must be the same")
    result = constructObjectTemplatedOnBasisAndResult(
        "DefaultIterativeSolver", basisFunctionType, resultType,
        linearOperator, gridFunction)
    result._space = linearOperator.domain()
    return result

%}

#endif
