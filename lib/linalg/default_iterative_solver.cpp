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

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "default_iterative_solver.hpp"

#include "real_wrapper_of_complex_thyra_linear_operator.hpp"
#include "../assembly/discrete_linear_operator.hpp"
#include "../assembly/linear_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <Teuchos_RCPBoostSharedPtrConversions.hpp>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::DefaultIterativeSolver(
        const LinearOperator<BasisFunctionType, ResultType>& linOp,
        const GridFunction<BasisFunctionType, ResultType>& rhsGridFun) :
    m_belosSolverWrapper(
        Teuchos::rcp<const Thyra::LinearOpBase<ResultType> >(
            linOp.weakForm())),
    m_space(linOp.domain()),
    // TODO: gridFun.projections should return a shared pointer to Vector
    // rather than a Vector
    m_rhs(new Vector<ResultType>(rhsGridFun.projections()))
{
    if (&linOp.dualToRange() != &rhsGridFun.dualSpace())
        throw std::runtime_error("DefaultIterativeSolver::DefaultIterativeSolver(): "
                                 "spaces do not match");
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::addPreconditioner(
        const Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> >& preconditioner)
{
    m_belosSolverWrapper.addPreconditioner(preconditioner);
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::initializeSolver(
        const Teuchos::RCP<Teuchos::ParameterList>& paramList)
{
    m_belosSolverWrapper.initializeSolver(paramList);
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::solve()
{
    assert(m_rhs->domain()->dim() == 1);

    const size_t size = m_rhs->range()->dim();
    m_sol.set_size(size);
    m_sol.fill(static_cast<ResultType>(0.));

    typedef Thyra::DefaultSpmdVector<ResultType> DenseVector;
    Teuchos::ArrayRCP<ResultType> solArray =
            Teuchos::arcp(m_sol.memptr(), 0 /* lowerOffset */,
                          size, false /* doesn't own memory */);
    DenseVector sol(Thyra::defaultSpmdVectorSpace<ResultType>(size),
                    solArray, 1 /* stride */);

    m_status = m_belosSolverWrapper.solve(
                Thyra::NOTRANS, *m_rhs, Teuchos::ptr(&sol));
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::getResult() const
{
    return gridFunctionFromCoefficients(m_space,
                                        m_space, // is this the right choice?
                                        m_sol);
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
Thyra::SolveStatus<typename DefaultIterativeSolver<BasisFunctionType, ResultType>::MagnitudeType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::getThyraSolveStatus() const
{
    return m_status;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);

} // namespace Bempp

#endif // WITH_TRILINOS
