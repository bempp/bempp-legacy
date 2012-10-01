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

#include "default_direct_solver.hpp"

#include "../assembly/abstract_boundary_operator.hpp"
#include "../assembly/blocked_boundary_operator.hpp"
#include "../assembly/boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/variant.hpp>

namespace Bempp
{

/** \cond HIDDEN_INTERNAL */

template <typename BasisFunctionType, typename ResultType> 
struct DefaultDirectSolver<BasisFunctionType, ResultType>::Impl
{
    Impl(const BoundaryOperator<BasisFunctionType, ResultType>& op_) :
        op(op_)
    {
    }

    Impl(const BlockedBoundaryOperator<BasisFunctionType, ResultType>& op_) :
        op(op_)
    {
    }
    
    boost::variant<
        BoundaryOperator<BasisFunctionType, ResultType>,
        BlockedBoundaryOperator<BasisFunctionType, ResultType> > op;
};

template <typename BasisFunctionType, typename ResultType>
DefaultDirectSolver<BasisFunctionType, ResultType>::DefaultDirectSolver(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp) :
    m_impl(new Impl(boundaryOp))
{
}

template <typename BasisFunctionType, typename ResultType>
DefaultDirectSolver<BasisFunctionType, ResultType>::DefaultDirectSolver(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp) :
    m_impl(new Impl(boundaryOp))
{
}

template <typename BasisFunctionType, typename ResultType>
DefaultDirectSolver<BasisFunctionType, ResultType>::~DefaultDirectSolver()
{
}

template <typename BasisFunctionType, typename ResultType>
Solution<BasisFunctionType, ResultType> 
DefaultDirectSolver<BasisFunctionType, ResultType>::solveImplNonblocked(
        const GridFunction<BasisFunctionType, ResultType>& rhs) const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;

    const BoundaryOp* boundaryOp = boost::get<BoundaryOp>(&m_impl->op);
    if (!boundaryOp)
        throw std::logic_error(
            "DefaultDirectSolver::solve(): for solvers constructed "
            "from a BlockedBoundaryOperator the other solve() overload "
            "must be used");
    Solver<BasisFunctionType, ResultType>::checkConsistency(
        *boundaryOp, rhs, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

    arma::Col<ResultType> armaSolution = arma::solve(
                boundaryOp->weakForm()->asMatrix(),
                rhs.projections());
    
    return Solution<BasisFunctionType, ResultType>(
        GridFunction<BasisFunctionType, ResultType>(
            boundaryOp->context(), 
            boundaryOp->domain(), boundaryOp->domain(), // is this the right choice?
            armaSolution, GridFunction<BasisFunctionType, ResultType>::COEFFICIENTS),
        SolutionStatus::CONVERGED,
        SolutionBase<BasisFunctionType, ResultType>::unknownTolerance(),
        "Solver finished");
}

template <typename BasisFunctionType, typename ResultType>
BlockedSolution<BasisFunctionType, ResultType>
DefaultDirectSolver<BasisFunctionType, ResultType>::solveImplBlocked(
    const std::vector<GridFunction<BasisFunctionType, ResultType> >& rhs) const
{
    typedef BlockedBoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;

    const BoundaryOp* boundaryOp = boost::get<BoundaryOp>(&m_impl->op);
    if (!boundaryOp)
        throw std::logic_error(
            "DefaultDirectSolver::solve(): for solvers constructed "
            "from a (non-blocked) BoundaryOperator the other solve() overload "
            "must be used");
    std::vector<GridFunction<BasisFunctionType, ResultType> > canonicalRhs =
            Solver<BasisFunctionType, ResultType>::canonicalizeBlockedRhs(
                *boundaryOp, rhs,
                ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
    // Shouldn't be needed, but better safe than sorry...
    Solver<BasisFunctionType, ResultType>::checkConsistency(
                *boundaryOp, canonicalRhs,
                ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

    // Construct the right-hand size vector
    arma::Col<ResultType> armaRhs(boundaryOp->totalGlobalDofCountInDualsToRanges());
    for (size_t i = 0, start = 0; i < canonicalRhs.size(); ++i) {
        const arma::Col<ResultType>& chunkProjections = canonicalRhs[i].projections();
        size_t chunkSize = chunkProjections.n_rows;
        armaRhs.rows(start, start + chunkSize - 1) = chunkProjections;
        start += chunkSize;
    }        

    // Solve
    arma::Col<ResultType> armaSolution = arma::solve(
                boundaryOp->weakForm()->asMatrix(),
                armaRhs);
    
    // Convert chunks of the solution vector into grid functions
    std::vector<GridFunction<BasisFunctionType, ResultType> > solutionFunctions;
    Solver<BasisFunctionType, ResultType>::constructBlockedGridFunction(
        armaSolution, *boundaryOp, solutionFunctions);

    // Return solution
    return BlockedSolution<BasisFunctionType, ResultType>(
        solutionFunctions,
        SolutionStatus::CONVERGED,
        SolutionBase<BasisFunctionType, ResultType>::unknownTolerance(),
        "Solver finished");
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultDirectSolver);

} // namespace Bempp
