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

#include "bempp/common/config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "default_iterative_solver.hpp"

#include "belos_solver_wrapper.hpp"
#include "solution.hpp"
#include "blocked_solution.hpp"
#include "../assembly/abstract_boundary_operator.hpp"
#include "../assembly/abstract_boundary_operator_pseudoinverse.hpp"
#include "../assembly/blocked_boundary_operator.hpp"
#include "../assembly/boundary_operator.hpp"
#include "../assembly/context.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator_composition.hpp"
#include "../assembly/identity_operator.hpp"
#include "../assembly/vector.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"

#include <Teuchos_RCPBoostSharedPtrConversions.hpp>
#include <Thyra_DefaultSpmdVectorSpace.hpp>

#include <boost/make_shared.hpp>
#include <boost/variant.hpp>

#include <tbb/task_scheduler_init.h>

namespace Bempp {

template <typename ValueType>
Teuchos::RCP<Thyra::DefaultSpmdVector<ValueType>>
wrapInTrilinosVector(arma::Col<ValueType> &col) {
  size_t size = col.n_rows;
  Teuchos::ArrayRCP<ValueType> trilinosArray = Teuchos::arcp(
      col.memptr(), 0 /* lowerOffset */, size, false /* doesn't own memory */);
  typedef Thyra::DefaultSpmdVector<ValueType> TrilinosVector;
  return Teuchos::RCP<TrilinosVector>(
      new TrilinosVector(Thyra::defaultSpmdVectorSpace<ValueType>(size),
                         trilinosArray, 1 /* stride */));
}

/** \cond HIDDEN_INTERNAL */

template <typename BasisFunctionType, typename ResultType>
struct DefaultIterativeSolver<BasisFunctionType, ResultType>::Impl {

  // Constructor for non-blocked operators
  Impl(const BoundaryOperator<BasisFunctionType, ResultType> &op_,
       ConvergenceTestMode::Mode mode_)
      : op(op_), mode(mode_) {
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    typedef Solver<BasisFunctionType, ResultType> Solver_;
    const BoundaryOp &boundaryOp = boost::get<BoundaryOp>(op);
    if (!boundaryOp.isInitialized())
      throw std::invalid_argument("DefaultIterativeSolver::Impl::Impl(): "
                                  "boundary operator must be initialized");

    if (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE) {
      if (boundaryOp.domain()->globalDofCount() !=
          boundaryOp.dualToRange()->globalDofCount())
        throw std::invalid_argument("DefaultIterativeSolver::Impl::Impl(): "
                                    "non-square system provided");

      solverWrapper.reset(new BelosSolverWrapper<ResultType>(
          Teuchos::rcp<const Thyra::LinearOpBase<ResultType>>(
              boundaryOp.weakForm())));
    } else if (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_RANGE) {
      if (boundaryOp.domain()->globalDofCount() !=
          boundaryOp.range()->globalDofCount())
        throw std::invalid_argument("DefaultIterativeSolver::Impl::Impl(): "
                                    "non-square system provided");

      BoundaryOp id =
          identityOperator(boundaryOp.context(), boundaryOp.range(),
                           boundaryOp.range(), boundaryOp.dualToRange());
      pinvId = pseudoinverse(id, boundaryOp.dualToRange());
      // dualToRange could be anything here.
      shared_ptr<DiscreteBoundaryOperator<ResultType>> totalBoundaryOp =
          boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType>>(
              boost::get<BoundaryOp>(pinvId).weakForm(), boundaryOp.weakForm());
      solverWrapper.reset(new BelosSolverWrapper<ResultType>(
          Teuchos::rcp<const Thyra::LinearOpBase<ResultType>>(
              totalBoundaryOp)));
    } else
      throw std::invalid_argument(
          "DefaultIterativeSolver::DefaultIterativeSolver(): "
          "invalid convergence test mode");
  }

  // Constructor for blocked operators
  Impl(const BlockedBoundaryOperator<BasisFunctionType, ResultType> &op_,
       ConvergenceTestMode::Mode mode_)
      : op(op_), mode(mode_) {
    typedef BlockedBoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    typedef Solver<BasisFunctionType, ResultType> Solver_;
    const BoundaryOp &boundaryOp = boost::get<BoundaryOp>(op);

    if (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE) {
      if (boundaryOp.totalGlobalDofCountInDomains() !=
          boundaryOp.totalGlobalDofCountInDualsToRanges())
        throw std::invalid_argument("DefaultIterativeSolver::Impl::Impl(): "
                                    "non-square system provided");
      solverWrapper.reset(new BelosSolverWrapper<ResultType>(
          Teuchos::rcp<const Thyra::LinearOpBase<ResultType>>(
              boundaryOp.weakForm())));
    } else if (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_RANGE) {
      if (boundaryOp.totalGlobalDofCountInDomains() !=
          boundaryOp.totalGlobalDofCountInRanges())
        throw std::invalid_argument("DefaultIterativeSolver::Impl::Impl(): "
                                    "non-square system provided");

      // Construct a block-diagonal operator composed of pseudoinverses
      // for appropriate spaces
      BlockedOperatorStructure<BasisFunctionType, ResultType> pinvIdStructure;
      size_t rowCount = boundaryOp.rowCount();
      size_t columnCount = boundaryOp.columnCount();
      for (size_t row = 0; row < rowCount; ++row) {
        // Find the first non-zero block in row #row and retrieve its context
        shared_ptr<const Context<BasisFunctionType, ResultType>> context;
        for (int col = 0; col < columnCount; ++col)
          if (boundaryOp.block(row, col).context()) {
            context = boundaryOp.block(row, col).context();
            break;
          }
        assert(context);

        BoundaryOperator<BasisFunctionType, ResultType> id = identityOperator(
            context, boundaryOp.range(row), boundaryOp.range(row),
            boundaryOp.dualToRange(row));
        pinvIdStructure.setBlock(row, row, pseudoinverse(id));
      }
      pinvId = BlockedBoundaryOperator<BasisFunctionType, ResultType>(
          pinvIdStructure);

      shared_ptr<DiscreteBoundaryOperator<ResultType>> totalBoundaryOp =
          boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType>>(
              boost::get<BoundaryOp>(pinvId).weakForm(), boundaryOp.weakForm());
      solverWrapper.reset(new BelosSolverWrapper<ResultType>(
          Teuchos::rcp<const Thyra::LinearOpBase<ResultType>>(
              totalBoundaryOp)));
    } else
      throw std::invalid_argument(
          "DefaultIterativeSolver::DefaultIterativeSolver(): "
          "invalid convergence test mode");
  }

  boost::variant<BoundaryOperator<BasisFunctionType, ResultType>,
                 BlockedBoundaryOperator<BasisFunctionType, ResultType>> op;
  ConvergenceTestMode::Mode mode;
  boost::scoped_ptr<BelosSolverWrapper<ResultType>> solverWrapper;
  boost::variant<BoundaryOperator<BasisFunctionType, ResultType>,
                 BlockedBoundaryOperator<BasisFunctionType, ResultType>> pinvId;
};

/** \endcond */

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::DefaultIterativeSolver(
    const BoundaryOperator<BasisFunctionType, ResultType> &boundaryOp,
    ConvergenceTestMode::Mode mode)
    : m_impl(new Impl(boundaryOp, mode)) {}

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::DefaultIterativeSolver(
    const BlockedBoundaryOperator<BasisFunctionType, ResultType> &boundaryOp,
    ConvergenceTestMode::Mode mode)
    : m_impl(new Impl(boundaryOp, mode)) {}

template <typename BasisFunctionType, typename ResultType>
DefaultIterativeSolver<BasisFunctionType,
                       ResultType>::~DefaultIterativeSolver() {}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::setPreconditioner(
    const Preconditioner<ResultType> &preconditioner) {
  m_impl->solverWrapper->setPreconditioner(preconditioner.get());
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::initializeSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &paramList) {
  m_impl->solverWrapper->initializeSolver(paramList);
}

template <typename BasisFunctionType, typename ResultType>
void DefaultIterativeSolver<BasisFunctionType, ResultType>::initializeSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &paramList,
    const Preconditioner<ResultType> &preconditioner) {
  m_impl->solverWrapper->setPreconditioner(preconditioner.get());
  m_impl->solverWrapper->initializeSolver(paramList);
}

template <typename BasisFunctionType, typename ResultType>
Solution<BasisFunctionType, ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::solveImplNonblocked(
    const GridFunction<BasisFunctionType, ResultType> &rhs) const {
  typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
  typedef typename ScalarTraits<ResultType>::RealType MagnitudeType;
  typedef Thyra::MultiVectorBase<ResultType> TrilinosVector;

  const BoundaryOp *boundaryOp = boost::get<BoundaryOp>(&m_impl->op);
  if (!boundaryOp)
    throw std::logic_error(
        "DefaultIterativeSolver::solve(): for solvers constructed "
        "from a BlockedBoundaryOperator the other solve() overload "
        "must be used");
  Solver<BasisFunctionType, ResultType>::checkConsistency(*boundaryOp, rhs,
                                                          m_impl->mode);

  // Construct rhs vector
  Vector<ResultType> projectionsVector(
      rhs.projections(boundaryOp->dualToRange()));
  Teuchos::RCP<TrilinosVector> rhsVector;
  if (m_impl->mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE)
    rhsVector = Teuchos::rcpFromRef(projectionsVector);
  else {
    const size_t size = boundaryOp->range()->globalDofCount();
    rhsVector.reset(new Vector<ResultType>(size));
    boost::get<BoundaryOp>(m_impl->pinvId).weakForm()->apply(
        Thyra::NOTRANS, projectionsVector, rhsVector.ptr(), 1., 0.);
  }

  // Construct solution vector
  arma::Col<ResultType> armaSolution(rhsVector->range()->dim());
  armaSolution.fill(static_cast<ResultType>(0.));
  Teuchos::RCP<TrilinosVector> solutionVector =
      wrapInTrilinosVector(armaSolution);

  // Get number of threads
  Fiber::ParallelizationOptions parallelOptions =
      boundaryOp->context()->assemblyOptions().parallelizationOptions();
  int maxThreadCount = 1;
  if (!parallelOptions.isOpenClEnabled()) {
    if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
      maxThreadCount = tbb::task_scheduler_init::automatic;
    else
      maxThreadCount = parallelOptions.maxThreadCount();
  }

  // Solve
  Thyra::SolveStatus<MagnitudeType> status;
  {
    // Initialize TBB threads here (to prevent their construction and
    // destruction on every matrix-vector multiplication)
    tbb::task_scheduler_init scheduler(maxThreadCount);
    status = m_impl->solverWrapper->solve(Thyra::NOTRANS, *rhsVector,
                                          solutionVector.ptr());
  }

  // Construct grid function and return
  return Solution<BasisFunctionType, ResultType>(
      GridFunction<BasisFunctionType, ResultType>(
          boundaryOp->context(), boundaryOp->domain(), armaSolution),
      status);
}

template <typename BasisFunctionType, typename ResultType>
BlockedSolution<BasisFunctionType, ResultType>
DefaultIterativeSolver<BasisFunctionType, ResultType>::solveImplBlocked(
    const std::vector<GridFunction<BasisFunctionType, ResultType>> &rhs) const {
  typedef BlockedBoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
  typedef typename ScalarTraits<ResultType>::RealType MagnitudeType;
  typedef Thyra::MultiVectorBase<ResultType> TrilinosVector;

  const BoundaryOp *boundaryOp = boost::get<BoundaryOp>(&m_impl->op);
  if (!boundaryOp)
    throw std::logic_error(
        "DefaultIterativeSolver::solve(): for solvers constructed "
        "from a (non-blocked) BoundaryOperator the other solve() overload "
        "must be used");
  std::vector<GridFunction<BasisFunctionType, ResultType>> canonicalRhs =
      Solver<BasisFunctionType, ResultType>::canonicalizeBlockedRhs(
          *boundaryOp, rhs, m_impl->mode);
  // Shouldn't be needed, but better safe than sorry...
  Solver<BasisFunctionType, ResultType>::checkConsistency(
      *boundaryOp, canonicalRhs, m_impl->mode);

  // Currently we only support convergence testing in space dual to range.

  // Construct the right-hand-side vector
  arma::Col<ResultType> armaProjections(
      boundaryOp->totalGlobalDofCountInDualsToRanges());
  for (size_t i = 0, start = 0; i < canonicalRhs.size(); ++i) {
    const arma::Col<ResultType> &chunkProjections =
        canonicalRhs[i].projections(boundaryOp->dualToRange(i));
    size_t chunkSize = chunkProjections.n_rows;
    armaProjections.rows(start, start + chunkSize - 1) = chunkProjections;
    start += chunkSize;
  }

  Vector<ResultType> projectionsVector(armaProjections);
  Teuchos::RCP<TrilinosVector> rhsVector;
  if (m_impl->mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE)
    rhsVector = Teuchos::rcpFromRef(projectionsVector);
  else {
    const size_t rhsSize = boundaryOp->totalGlobalDofCountInRanges();
    rhsVector.reset(new Vector<ResultType>(rhsSize));
    boost::get<BoundaryOp>(m_impl->pinvId).weakForm()->apply(
        Thyra::NOTRANS, projectionsVector, rhsVector.ptr(), 1., 0.);
  }

  // Initialize the solution vector
  size_t solutionSize = 0;
  for (size_t i = 0; i < canonicalRhs.size(); ++i)
    solutionSize += boundaryOp->domain(i)->globalDofCount();
  arma::Col<ResultType> armaSolution(solutionSize);
  armaSolution.fill(static_cast<ResultType>(0.));
  Teuchos::RCP<TrilinosVector> solutionVector =
      wrapInTrilinosVector(armaSolution);

  // Get context of the first non-empty operator
  size_t rowCount = boundaryOp->rowCount();
  shared_ptr<const Context<BasisFunctionType, ResultType>> context;
  for (size_t row = 0; row < rowCount; ++row)
    if (boundaryOp->block(row, 0).context()) {
      context = boundaryOp->block(row, 0).context();
      break;
    }
  assert(context);

  // Get number of threads
  Fiber::ParallelizationOptions parallelOptions =
      context->assemblyOptions().parallelizationOptions();
  int maxThreadCount = 1;
  if (!parallelOptions.isOpenClEnabled()) {
    if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
      maxThreadCount = tbb::task_scheduler_init::automatic;
    else
      maxThreadCount = parallelOptions.maxThreadCount();
  }

  // Solve
  Thyra::SolveStatus<MagnitudeType> status;
  {
    // Initialize TBB threads here (to prevent their construction and
    // destruction on every matrix-vector multiplication)
    tbb::task_scheduler_init scheduler(maxThreadCount);
    status = m_impl->solverWrapper->solve(Thyra::NOTRANS, *rhsVector,
                                          solutionVector.ptr());
  }

  // Convert chunks of the solution vector into grid functions
  std::vector<GridFunction<BasisFunctionType, ResultType>> solutionFunctions;
  Solver<BasisFunctionType, ResultType>::constructBlockedGridFunction(
      armaSolution, *boundaryOp, solutionFunctions);

  // Return solution
  return BlockedSolution<BasisFunctionType, ResultType>(solutionFunctions,
                                                        status);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultIterativeSolver);

} // namespace Bempp

#endif // WITH_TRILINOS
