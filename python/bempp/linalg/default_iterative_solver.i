%{
#include "linalg/default_iterative_solver.hpp"
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
    return constructObjectTemplatedOnBasisAndResult(
        "DefaultIterativeSolver", basisFunctionType, resultType,
        linearOperator, gridFunction)

%}

#endif
