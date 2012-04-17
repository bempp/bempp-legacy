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

template<typename ValueType>
DefaultIterativeSolver<ValueType>::DefaultIterativeSolver(
        const LinearOperator<ValueType>& linOp,
        const GridFunction<ValueType>& gridFun) :
    m_discreteOperator(linOp.assembledDiscreteLinearOperator()),
    m_rhs(new Vector<ValueType>(gridFun.coefficients())),
    m_space(linOp.trialSpace())
{
    if (!linOp.isAssembled())
        throw std::runtime_error("DefaultIterativeSolver::DefaultIterativeSolver(): "
                                 "operator is not assembled");
    if (&linOp.trialSpace() != &gridFun.space())
        throw std::runtime_error("DefaultIterativeSolver::DefaultIterativeSolver(): "
                                 "spaces do not match");
}

template<typename ValueType>
void DefaultIterativeSolver<ValueType>::addPreconditioner(
        Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > preconditioner)
{
    m_preconditioner = preconditioner;
}

template<typename ValueType>
void DefaultIterativeSolver<ValueType>::initializeSolver(Teuchos::RCP<Teuchos::ParameterList> paramList)
{
    Teuchos::RCP<Teuchos::FancyOStream> out =
            Teuchos::VerboseObjectBase::getDefaultOStream();

    Thyra::BelosLinearOpWithSolveFactory<ValueType> invertibleOpFactory;
    invertibleOpFactory.setParameterList(paramList);
    invertibleOpFactory.setOStream(out);
    invertibleOpFactory.setVerbLevel(Teuchos::VERB_DEFAULT);
    // Wrap discreteOperator in reference counted Trilinos pointer.
    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > trilinosDiscreteLhs(
                &m_discreteOperator, false /*don't own*/);
    Teuchos::RCP<const Thyra::LinearOpSourceBase<ValueType> > linearOpSourcePtr(
                new Thyra::DefaultLinearOpSource<ValueType>(trilinosDiscreteLhs));

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

template <typename ValueType>
void DefaultIterativeSolver<ValueType>::solve()
{
    const size_t size = m_rhs->range()->dim();
    const size_t nrhs = m_rhs->domain()->dim();

    if (m_lhs.is_null()) throw std::runtime_error("Solver not initialized");

    m_sol = Teuchos::RCP<Thyra::MultiVectorBase<ValueType> >(new Thyra::DefaultSpmdMultiVector<ValueType>(
                Thyra::defaultSpmdVectorSpace<ValueType>(size),
                Thyra::defaultSpmdVectorSpace<ValueType>(nrhs)));
    Thyra::assign(m_sol.ptr(), 0.);

    m_status = m_lhs->solve(Thyra::NOTRANS, *m_rhs, m_sol.ptr());
}

template <typename ValueType>
GridFunction<ValueType> DefaultIterativeSolver<ValueType>::getResult() const
{
    const size_t size = m_sol->range()->dim();
    const size_t nrhs = m_sol->domain()->dim();

    Thyra::ConstDetachedMultiVectorView<ValueType> solView(m_sol);
    return GridFunction<ValueType>(
                m_space, arma::Mat<ValueType>(solView.values(), size, nrhs));
}


template <typename ValueType>
typename Solver<ValueType>::EStatus DefaultIterativeSolver<ValueType>::getStatus() const
{
    switch (m_status.solveStatus) {
    case Thyra::SOLVE_STATUS_CONVERGED:
        return Solver<ValueType>::CONVERGED;
    case Thyra::SOLVE_STATUS_UNCONVERGED:
        return Solver<ValueType>::UNCONVERGED;
    default:
        return Solver<ValueType>::UNKNOWN;
    }
}

template <typename ValueType>
double DefaultIterativeSolver<ValueType>::getSolveTolerance() const
{
    return m_status.achievedTol;
}

template <typename ValueType>
std::string DefaultIterativeSolver<ValueType>::getSolverMessage() const
{
    return m_status.message;
}

template <typename ValueType>
Thyra::SolveStatus<ValueType>
DefaultIterativeSolver<ValueType>::getThyraSolveStatus() const
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

#ifdef COMPILE_FOR_FLOAT
template class DefaultIterativeSolver<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class DefaultIterativeSolver<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class DefaultIterativeSolver<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class DefaultIterativeSolver<std::complex<double> >;
#endif

} // namespace Bempp

#endif // WITH_TRILINOS
