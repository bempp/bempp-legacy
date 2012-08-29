%{
#include "linalg/default_iterative_solver.hpp"
#include "linalg/blocked_solution.hpp" // temp
#include "assembly/blocked_boundary_operator.hpp" // temp
#include "linalg/belos_solver_wrapper.hpp"
#include "linalg/solver.hpp"
#include "bempp/common/config_trilinos.hpp"
#include "assembly/boundary_operator.hpp"
#include "linalg/solver.hpp"


%}

// TODO
// %include "default_iterative_solver_docstrings.i"

#ifdef WITH_TRILINOS




namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(
    double tol, int maxIterationCount = 1000);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(
    double tol, int maxIterationCount = 1000);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "linalg/default_iterative_solver.hpp"
#undef shared_ptr

//namespace Bempp {

//    template <typename BasisFunctionType, typename ResultType>
//    class DefaultIterativeSolver : public Solver<BasisFunctionType, ResultType>
//    {
//    public:

//        DefaultIterativeSolver(
//                const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
//                ConvergenceTestMode::Mode mode=ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

//        DefaultIterativeSolver(
//                const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
//                ConvergenceTestMode::Mode mode=ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

//        virtual ~DefaultIterativeSolver();

//        void setPreconditioner(
//                const Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> >& preconditioner);

//        void initializeSolver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

//    };



//}






namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);
} // namespace Bempp

%pythoncode %{

def defaultIterativeSolver(boundaryOperator,test_convergence):
    """Construct the default linear solver.

    This solver lets you solve the equation A f = g for the function
    f, with A being the boundary operator passed via the boundaryOperator
    argument and g a grid function supplied to the solve() method.
    """
    basisFunctionType = boundaryOperator.basisFunctionType()
    resultType = boundaryOperator.resultType()
    result = constructObjectTemplatedOnBasisAndResult(
        "DefaultIterativeSolver", basisFunctionType, resultType,
        boundaryOperator,test_convergence)
    result._boundaryOperator = boundaryOperator
    return result

%}

#endif
