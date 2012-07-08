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

#define shared_ptr boost::shared_ptr
%include "linalg/default_iterative_solver.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);
} // namespace Bempp

%pythoncode %{

def defaultIterativeSolver(boundaryOperator, gridFunction):
    """Construct the default linear solver.

    This solver lets you solve the equation A f = g for the function
    f, with A being the boundary operator passed via the boundaryOperator
    argument and g the function represented by the gridFunction
    argument.
    """
    basisFunctionType = boundaryOperator.basisFunctionType()
    if (basisFunctionType != gridFunction.basisFunctionType()):
        raise TypeError("BasisFunctionType of boundaryOperator and "
                        "gridFunction must be the same")
    resultType = boundaryOperator.resultType()
    if (resultType != gridFunction.resultType()):
        raise TypeError("ResultType of boundaryOperator and "
                        "gridFunction must be the same")
    result = constructObjectTemplatedOnBasisAndResult(
        "DefaultIterativeSolver", basisFunctionType, resultType,
        boundaryOperator, gridFunction)
    result._space = boundaryOperator.domain()
    return result

%}

#endif
