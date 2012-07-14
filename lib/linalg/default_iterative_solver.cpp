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

#include "belos_solver_wrapper.hpp"
#include "solution.hpp"
#include "blocked_solution.hpp"
#include "../assembly/abstract_boundary_operator.hpp"
#include "../assembly/abstract_boundary_operator_pseudoinverse.hpp"
#include "../assembly/blocked_boundary_operator.hpp"
#include "../assembly/boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator_composition.hpp"
#include "../assembly/identity_operator.hpp"
#include "../assembly/vector.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"

#include <Teuchos_RCPBoostSharedPtrConversions.hpp>
#include <Thyra_DefaultSpmdVectorSpace.hpp>

#include <boost/make_shared.hpp>
#include <boost/variant.hpp>

namespace Bempp
{

template <typename ValueType>
Teuchos::RCP<Thyra::DefaultSpmdVector<ValueType> >
wrapInTrilinosVector(arma::Col<ValueType>& col)
{
    size_t size = col.n_rows;
    Teuchos::ArrayRCP<ValueType> trilinosArray =
            Teuchos::arcp(col.memptr(), 0 /* lowerOffset */,
                          size, false /* doesn't own memory */);
    typedef Thyra::DefaultSpmdVector<ValueType> TrilinosVector;
    return Teuchos::RCP<TrilinosVector>(new TrilinosVector(
        Thyra::defaultSpmdVectorSpace<ValueType>(size),
        trilinosArray, 1 /* stride */));
}

template <typename BasisFunctionType, typename ResultType> 
struct DefaultIterativeSolver<BasisFunctionType, ResultType>::Impl
{
    Impl(const BoundaryOperator<BasisFunctionType, ResultType>& op_,
         ConvergenceTestMode mode_) :
        op(op_),
        mode(mode_)
    {
        typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
        const BoundaryOp& boundaryOp = boost::get<BoundaryOp>(op);

        if (mode == TEST_CONVERGENCE_IN_DUAL_TO_RANGE) {
            solverWrapper.reset(
                        new BelosSolverWrapper<ResultType>(
                            Teuchos::rcp<const Thyra::LinearOpBase<ResultType> >(
                                boundaryOp.weakForm())));
            // m_rhs.reset(new Vector<ResultType>(rhsGridFun.projections()));
        }
        else if (mode == TEST_CONVERGENCE_IN_RANGE) {
            BoundaryOp id = identityOperator(
                        boundaryOp.context(), boundaryOp.range(), boundaryOp.range(),
                        boundaryOp.dualToRange(), "I");
            pinvId = pseudoinverse(id);
            shared_ptr<DiscreteBoundaryOperator<ResultType> > totalBoundaryOp =
                    boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                        pinvId.weakForm(),
                        boundaryOp.weakForm());
            solverWrapper.reset(
                        new BelosSolverWrapper<ResultType>(
                            Teuchos::rcp<const Thyra::LinearOpBase<ResultType> >(
                                totalBoundaryOp)));
        }
        else
            throw std::invalid_argument(
                    "DefaultIterativeSolver::DefaultIterativeSolver(): "
                    "invalid convergence test mode");
    }

    Impl(const BlockedBoundaryOperator<BasisFunctionType, ResultType>& op_,
         ConvergenceTestMode mode_) :
        op(op_),
        mode(mode_),
        solverWrapper(new BelosSolverWrapper<ResultType>(
                          Teuchos::rcp<const Thyra::LinearOpBase<ResultType> >(
                              op_.weakForm())))
    {
        if (mode != TEST_CONVERGENCE_IN_DUAL_TO_RANGE)
            throw std::runtime_error(
                "DefaultIterativeSolver::DefaultIterativeSolver(): "
                "currently equations involving blocked operators can only "
                "be solved with convergence tested in space dual to range");
    }
    
    boost::variant<
        BoundaryOperator<BasisFunctionType, ResultType>,
        BlockedBoundaryOperator<BasisFunctionType, ResultType> > op;
    typename DefaultIterativeSolver<BasisFunctionType, ResultType>::
    ConvergenceTestMode mode;
    boost::scoped_ptr<BelosSolverWrapper<ResultType> > solverWrapper;
    BoundaryOperator<BasisFunctionType, ResultType> pinvId;
};

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::DefaultIterativeSolver(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        ConvergenceTestMode mode) :
    m_impl(new Impl(boundaryOp, mode))
{
}

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::DefaultIterativeSolver(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        ConvergenceTestMode mode) :
    m_impl(new Impl(boundaryOp, mode))
{
}

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::~DefaultIterativeSolver()
{
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::setPreconditioner(
        const Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> >& preconditioner)
{
    m_impl->solverWrapper->setPreconditioner(preconditioner);
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::initializeSolver(
        const Teuchos::RCP<Teuchos::ParameterList>& paramList)
{
    m_impl->solverWrapper->initializeSolver(paramList);
}

template <typename BasisFunctionType, typename ResultType>
Solution<BasisFunctionType, ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::solve(
        const GridFunction<BasisFunctionType, ResultType>& rhs) const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    typedef Thyra::MultiVectorBase<ResultType> TrilinosVector;

    const BoundaryOp* boundaryOp = boost::get<BoundaryOp>(&m_impl->op);
    if (!boundaryOp)
        throw std::logic_error(
            "DefaultIterativeSolver::solve(): for solvers constructed "
            "from a BlockedBoundaryOperator the other solve() overload "
            "must be used");
    if ((m_impl->mode == TEST_CONVERGENCE_IN_DUAL_TO_RANGE &&
         rhs.dualSpace() != boundaryOp->dualToRange()) ||
        (m_impl->mode == TEST_CONVERGENCE_IN_RANGE &&
         rhs.space() != boundaryOp->range()))
        throw std::invalid_argument(
            "DefaultIterativeSolver::solve(): spaces do not match");

    // Construct rhs vector
    Vector<ResultType> projectionsVector(rhs.projections());
    Teuchos::RCP<TrilinosVector> rhsVector;
    if (m_impl->mode == TEST_CONVERGENCE_IN_DUAL_TO_RANGE)
        rhsVector = Teuchos::rcpFromRef(projectionsVector);
    else {
        const size_t size = boundaryOp->range()->globalDofCount();
        rhsVector.reset(new Vector<ResultType>(size));
        m_impl->pinvId.weakForm()->apply(
            Thyra::NOTRANS, projectionsVector, rhsVector.ptr(), 1., 0.);
    }

    // Construct solution vector
    arma::Col<ResultType> armaSolution(rhsVector->range()->dim());
    armaSolution.fill(static_cast<ResultType>(0.));
    Teuchos::RCP<TrilinosVector> solutionVector = wrapInTrilinosVector(armaSolution);

    // Solve
    Thyra::SolveStatus<MagnitudeType> status = m_impl->solverWrapper->solve(
        Thyra::NOTRANS, *rhsVector, solutionVector.ptr());

    // Construct grid function and return
    return Solution<BasisFunctionType, ResultType>(
        GridFunction<BasisFunctionType, ResultType>(
            boundaryOp->context(), 
            boundaryOp->domain(), boundaryOp->domain(), // is this the right choice?
            armaSolution, GridFunction<BasisFunctionType, ResultType>::COEFFICIENTS),
        status);
}

template <typename BasisFunctionType, typename ResultType>
BlockedSolution<BasisFunctionType, ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::solve(
    const std::vector<GridFunction<BasisFunctionType, ResultType> >& rhs) const
{
    typedef BlockedBoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    const BoundaryOp* boundaryOp = boost::get<BoundaryOp>(&m_impl->op);

    if (!boundaryOp)
        throw std::logic_error(
            "DefaultIterativeSolver::solve(): for solvers constructed "
            "from a (non-blocked) BoundaryOperator the other solve() overload "
            "must be used");

    const size_t columnCount = boundaryOp->columnCount();
    const size_t rowCount = boundaryOp->rowCount();

    if (rhs.size() != rowCount)
        throw std::invalid_argument(
            "DefaultIterativeSolver::solve(): incorrect number of grid functions");
    for (size_t i = 0; i < rhs.size(); ++i)
        if (rhs[i].dualSpace() != boundaryOp->dualToRange(i))
            throw std::invalid_argument(
                "DefaultIterativeSolver::solve(): dual space of grid function #" +
                toString(i) + 
                " does not match that of the blocked boundary operator");

    typedef Thyra::DefaultSpmdVector<ResultType> DenseVector;

    // Currently we only support convergence testing in space dual to range.

    // Construct the right-hand-side vector
    arma::Col<ResultType> armaRhs(boundaryOp->totalGlobalDofCountInDualsToRanges());
    for (size_t i = 0, start = 0; i < rhs.size(); ++i) {
        const arma::Col<ResultType>& chunkProjections = rhs[i].projections();
        size_t chunkSize = chunkProjections.n_rows;
        armaRhs.rows(start, start + chunkSize - 1) = chunkProjections;
        start += chunkSize;
    }        
    Teuchos::RCP<DenseVector> rhsVector = wrapInTrilinosVector(armaRhs);

    // Initialize the solution vector
    size_t solutionSize = 0;
    for (size_t i = 0; i < columnCount; ++i)
        solutionSize += boundaryOp->domain(i)->globalDofCount();
    arma::Col<ResultType> armaSolution(solutionSize);
    armaSolution.fill(static_cast<ResultType>(0.));
    Teuchos::RCP<DenseVector> solutionVector = wrapInTrilinosVector(armaSolution);

    Thyra::SolveStatus<MagnitudeType> status = m_impl->solverWrapper->solve(
        Thyra::NOTRANS, *rhsVector, solutionVector.ptr());
    
    // Convert chunks of the solution vector into grid functions
    std::vector<GridFunction<BasisFunctionType, ResultType> > solutionFunctions(
        columnCount);
    for (size_t i = 0, start = 0; i < solutionFunctions.size(); ++i) {        
        // Find the first non-zero block in column i and retrieve its context
        shared_ptr<const Context<BasisFunctionType, ResultType> > context;
        for (int row = 0; row < rowCount; ++row)
            if (boundaryOp->block(row, i).context()) {
                context = boundaryOp->block(row, i).context();
                break;
            }
        assert(context);

        size_t chunkSize = boundaryOp->domain(i)->globalDofCount();
        solutionFunctions[i] = 
        GridFunction<BasisFunctionType, ResultType>(
            context, 
            boundaryOp->domain(i), boundaryOp->domain(i), // is this the right choice?
            armaSolution.rows(start, start + chunkSize - 1),
            GridFunction<BasisFunctionType, ResultType>::COEFFICIENTS);
        start += chunkSize;
    }
    // Return solution
    return BlockedSolution<BasisFunctionType, ResultType>(solutionFunctions, status);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);

} // namespace Bempp

#endif // WITH_TRILINOS
