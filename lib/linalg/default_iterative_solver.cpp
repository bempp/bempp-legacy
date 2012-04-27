// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "default_iterative_solver.hpp"

#include "config_trilinos.hpp"

#include "../assembly/linear_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#ifdef WITH_TRILINOS
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_BelosLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::DefaultIterativeSolver(
        const LinearOperator<BasisFunctionType, ResultType>& linOp,
        const GridFunction<BasisFunctionType, ResultType>& gridFun) :
    m_discreteOperator(linOp.assembledDiscreteLinearOperator()),
    m_rhs(new Vector<ResultType>(gridFun.coefficients())),
    m_space(linOp.trialSpace())
{
    if (!linOp.isAssembled())
        throw std::runtime_error("DefaultIterativeSolver::DefaultIterativeSolver(): "
                                 "operator is not assembled");
    if (&linOp.trialSpace() != &gridFun.space())
        throw std::runtime_error("DefaultIterativeSolver::DefaultIterativeSolver(): "
                                 "spaces do not match");
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::addPreconditioner(
        Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> > preconditioner)
{
    m_preconditioner = preconditioner;
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::initializeSolver(
        Teuchos::RCP<Teuchos::ParameterList> paramList)
{
    Teuchos::RCP<Teuchos::FancyOStream> out =
            Teuchos::VerboseObjectBase::getDefaultOStream();

    Thyra::BelosLinearOpWithSolveFactory<ResultType> invertibleOpFactory;
    invertibleOpFactory.setParameterList(paramList);
    invertibleOpFactory.setOStream(out);
    invertibleOpFactory.setVerbLevel(Teuchos::VERB_DEFAULT);
    // Wrap discreteOperator in reference counted Trilinos pointer.
    Teuchos::RCP<const Thyra::LinearOpBase<ResultType> > trilinosDiscreteLhs(
                &m_discreteOperator, false /*don't own*/);
    Teuchos::RCP<const Thyra::LinearOpSourceBase<ResultType> > linearOpSourcePtr(
                new Thyra::DefaultLinearOpSource<ResultType>(trilinosDiscreteLhs));

    if (!m_preconditioner.is_null()) {
        // Preconditioner defined
        m_lhs = invertibleOpFactory.createOp();
        invertibleOpFactory.initializePreconditionedOp(linearOpSourcePtr,
                m_preconditioner,
                m_lhs.get(),
                Thyra::SUPPORT_SOLVE_UNSPECIFIED);
    } else {
        // No preconditioner
        m_lhs = Thyra::linearOpWithSolve(invertibleOpFactory, trilinosDiscreteLhs);
    }
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::solve()
{
    const size_t size = m_rhs->range()->dim();
    const size_t nrhs = m_rhs->domain()->dim();

    if (m_lhs.is_null()) throw std::runtime_error("Solver not initialized");

    m_sol = Teuchos::RCP<Thyra::MultiVectorBase<ResultType> >(
                new Thyra::DefaultSpmdMultiVector<ResultType>(
                Thyra::defaultSpmdVectorSpace<ResultType>(size),
                Thyra::defaultSpmdVectorSpace<ResultType>(nrhs)));
    Thyra::assign(m_sol.ptr(), static_cast<ResultType>(0.));

    m_status = m_lhs->solve(Thyra::NOTRANS, *m_rhs, m_sol.ptr());
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::getResult() const
{
    const size_t size = m_sol->range()->dim();
    const size_t nrhs = m_sol->domain()->dim();

    Thyra::ConstDetachedMultiVectorView<ResultType> solView(m_sol);
    return GridFunction<BasisFunctionType, ResultType>(
                m_space, arma::Mat<ResultType>(solView.values(), size, nrhs));
}

template <typename BasisFunctionType, typename ResultType>
typename Solver<BasisFunctionType, ResultType>::EStatus
DefaultIterativeSolver<BasisFunctionType, ResultType>::getStatus() const
{
    switch (m_status.solveStatus) {
    case Thyra::SOLVE_STATUS_CONVERGED:
        return Solver<BasisFunctionType, ResultType>::CONVERGED;
    case Thyra::SOLVE_STATUS_UNCONVERGED:
        return Solver<BasisFunctionType, ResultType>::UNCONVERGED;
    default:
        return Solver<BasisFunctionType, ResultType>::UNKNOWN;
    }
}

template <typename BasisFunctionType, typename ResultType>
typename DefaultIterativeSolver<BasisFunctionType, ResultType>::MagnitudeType
DefaultIterativeSolver<BasisFunctionType, ResultType>::getSolveTolerance() const
{
    return m_status.achievedTol;
}

template <typename BasisFunctionType, typename ResultType>
std::string DefaultIterativeSolver<BasisFunctionType, ResultType>::getSolverMessage() const
{
    return m_status.message;
}

template <typename BasisFunctionType, typename ResultType>
Thyra::SolveStatus<ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::getThyraSolveStatus() const
{
    return m_status;
}

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol)
{
    Teuchos::RCP<Teuchos::ParameterList> paramList(
                new Teuchos::ParameterList("DefaultParameters"));
    paramList->set("Solver Type", "Pseudo Block GMRES");
    Teuchos::ParameterList& solverTypesList = paramList->sublist("Solver Types");
    Teuchos::ParameterList& pseudoBlockGmresList =
            solverTypesList.sublist("Pseudo Block GMRES");
    pseudoBlockGmresList.set("Convergence Tolerance", tol);
    return paramList;
}

Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol)
{
    Teuchos::RCP<Teuchos::ParameterList> paramList(
                new Teuchos::ParameterList("DefaultParameters"));
    paramList->set("Solver Type", "Pseudo Block CG");
    Teuchos::ParameterList& solverTypesList = paramList->sublist("Solver Types");
    Teuchos::ParameterList& pseudoBlockCgList =
            solverTypesList.sublist("Pseudo Block CG");
    pseudoBlockCgList.set("Convergence Tolerance", tol);
    return paramList;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT_REAL_ONLY(DefaultIterativeSolver);

} // namespace Bempp

#endif // WITH_TRILINOS
