// Copyright (C) 2011 by the BEM++ Authors
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

#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"

#include "assembly/blocked_boundary_operator.hpp"
#include "assembly/blocked_operator_structure.hpp"
#include "assembly/context.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"
#include "common/scalar_traits.hpp"
#include "grid/grid_factory.hpp"
#include "linalg/default_iterative_solver.hpp"
#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_complex.hpp>

#include <armadillo>

using namespace Bempp;

// Tests

template <typename ValueType_>
class UnitFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = 1.;
    }
};

BOOST_AUTO_TEST_SUITE(DefaultIterativeSolver)

BOOST_AUTO_TEST_CASE_TEMPLATE(both_convergence_testing_strategies_agree_for_dual_equal_to_domain,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    std::auto_ptr<Grid> grid = GridFactory::importGmshGrid(
                params, "../examples/meshes/sphere-152.msh");

    PiecewiseLinearContinuousScalarSpace<BFT> space(*grid);
    space.assignDofs();
    shared_ptr<Space<BFT> > spacePtr = make_shared_from_ref(space);

    AssemblyOptions assemblyOptions;
    NumericalQuadratureStrategy<BFT, RT> quadStrategy;
    Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);
    shared_ptr<Context<BFT, RT> > contextPtr = make_shared_from_ref(context);

    BoundaryOperator<BFT, RT> slpOp = laplace3dSingleLayerBoundaryOperator<BFT, RT>(
                contextPtr, spacePtr, spacePtr, spacePtr);
    BoundaryOperator<BFT, RT> dlpOp = laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
                contextPtr, spacePtr, spacePtr, spacePtr);
    BoundaryOperator<BFT, RT> id = identityOperator<BFT, RT>(
                contextPtr, spacePtr, spacePtr, spacePtr);

    BoundaryOperator<BFT, RT> lhsOp = slpOp;
    BoundaryOperator<BFT, RT> rhsOp = -0.5 * id + dlpOp;

    GridFunction<BFT, RT> u(contextPtr, spacePtr, spacePtr,
                surfaceNormalIndependentFunction(UnitFunctor<RT>()));
    GridFunction<BFT, RT> rhs = rhsOp * u;

    typedef Bempp::DefaultIterativeSolver<BFT, RT> IterSolver;
    const RealType solverTol = 1e-6;

    IterSolver solverDual(
        lhsOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
    solverDual.initializeSolver(defaultGmresParameterList(solverTol));
    Solution<BFT, RT> solutionDual = solverDual.solve(rhs);
    arma::Col<RT> solutionVectorDual = solutionDual.gridFunction().coefficients();

    IterSolver solverPrimal(
        lhsOp, IterSolver::TEST_CONVERGENCE_IN_RANGE);
    solverPrimal.initializeSolver(defaultGmresParameterList(solverTol));
    Solution<BFT, RT> solutionPrimal = solverPrimal.solve(rhs);
    arma::Col<RT> solutionVectorPrimal = solutionPrimal.gridFunction().coefficients();

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    solutionVectorDual, solutionVectorPrimal, solverTol * 1000.));
}

template <typename BFT, typename RT>
class LaplaceDirichletFixture
{
public:
    LaplaceDirichletFixture() {
        GridParameters params;
        params.topology = GridParameters::TRIANGULAR;
        grid = GridFactory::importGmshGrid(
            params, "../examples/meshes/sphere-152.msh");

        shared_ptr<Space<BFT> > pwiseConstants(
            new Bempp::PiecewiseConstantScalarSpace<BFT>(*grid));
        shared_ptr<Space<BFT> > pwiseLinears(
            new PiecewiseLinearContinuousScalarSpace<BFT>(*grid));
        pwiseConstants->assignDofs();
        pwiseLinears->assignDofs();

        AssemblyOptions assemblyOptions;
        shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy( 
            new NumericalQuadratureStrategy<BFT, RT>);
        shared_ptr<Context<BFT, RT> > context(
            new Context<BFT, RT>(quadStrategy, assemblyOptions));
        
        BoundaryOperator<BFT, RT> slpOp = laplace3dSingleLayerBoundaryOperator<BFT, RT>(
            context, pwiseConstants, pwiseLinears, pwiseConstants);
        BoundaryOperator<BFT, RT> dlpOp = laplace3dDoubleLayerBoundaryOperator<BFT, RT>(
            context, pwiseLinears, pwiseLinears, pwiseConstants);
        BoundaryOperator<BFT, RT> id = identityOperator<BFT, RT>(
            context, pwiseLinears, pwiseLinears, pwiseConstants);
        
        lhsOp = slpOp;
        BoundaryOperator<BFT, RT> rhsOp = -0.5 * id + dlpOp;
        
        GridFunction<BFT, RT> u(context, pwiseLinears, pwiseLinears,
                                surfaceNormalIndependentFunction(UnitFunctor<RT>()));
        rhs = rhsOp * u;
    }

    std::auto_ptr<Grid> grid;
    BoundaryOperator<BFT, RT> lhsOp;
    GridFunction<BFT, RT> rhs;
};

BOOST_AUTO_TEST_CASE_TEMPLATE(boundary_operator_agrees_with_trivial_1x1_blocked_boundary_operator_for_convergence_testing_in_dual,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    typedef Bempp::DefaultIterativeSolver<BFT, RT> IterSolver;
    const RealType solverTol = 1e-5;

    LaplaceDirichletFixture<BFT, RT> fixture;

    arma::Col<RT> solutionVectorNonblocked;

    // Solve using nonblocked operator
    {
        IterSolver solver(
            fixture.lhsOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);       
        solver.initializeSolver(defaultGmresParameterList(solverTol));
        Solution<BFT, RT> solution = solver.solve(fixture.rhs);
        solutionVectorNonblocked = solution.gridFunction().coefficients();
    }

    // Solve using trivial (1 x 1) blocked operator
    {
        BlockedOperatorStructure<BFT, RT> structure;
        structure.setBlock(0, 0, fixture.lhsOp);
        
        BlockedBoundaryOperator<BFT, RT> lhsBlockedOp(structure);
        std::vector<GridFunction<BFT, RT> > blockedRhs(1);
        blockedRhs[0] = fixture.rhs;

        IterSolver solver(
            lhsBlockedOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);       
        solver.initializeSolver(defaultGmresParameterList(solverTol));
        BlockedSolution<BFT, RT> solution = solver.solve(blockedRhs);
        arma::Col<RT> solutionVectorBlocked = solution.gridFunction(0).coefficients();

        BOOST_CHECK(check_arrays_are_close<ValueType>(
                        solutionVectorNonblocked, solutionVectorBlocked, solverTol * 10));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(boundary_operator_agrees_with_diagonal_2x2_blocked_boundary_operator_for_convergence_testing_in_dual,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    typedef Bempp::DefaultIterativeSolver<BFT, RT> IterSolver;
    const RealType solverTol = 1e-5;

    LaplaceDirichletFixture<BFT, RT> fixture;

    arma::Col<RT> solutionVectorNonblocked;

    // Solve using nonblocked operator
    {
        IterSolver solver(
            fixture.lhsOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);       
        solver.initializeSolver(defaultGmresParameterList(solverTol));
        Solution<BFT, RT> solution = solver.solve(fixture.rhs);
        solutionVectorNonblocked = solution.gridFunction().coefficients();
    }

    // Solve using diagonal 2x2 ([A, 0; 0, A]) blocked operator
    {
        BlockedOperatorStructure<BFT, RT> structure;
        structure.setBlock(0, 0, fixture.lhsOp);
        structure.setBlock(1, 1, fixture.lhsOp);
        
        BlockedBoundaryOperator<BFT, RT> lhsBlockedOp(structure);
        std::vector<GridFunction<BFT, RT> > blockedRhs(2);
        blockedRhs[0] = fixture.rhs;
        blockedRhs[1] = 2. * fixture.rhs;

        IterSolver solver(
            lhsBlockedOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);       
        solver.initializeSolver(defaultGmresParameterList(solverTol));
        BlockedSolution<BFT, RT> solution = solver.solve(blockedRhs);
        arma::Col<RT> solutionVectorBlock0 = solution.gridFunction(0).coefficients();
        arma::Col<RT> solutionVectorBlock1 = solution.gridFunction(1).coefficients() / 2.;

        BOOST_CHECK(check_arrays_are_close<ValueType>(
                        solutionVectorNonblocked, solutionVectorBlock0, solverTol * 10));
        BOOST_CHECK(check_arrays_are_close<ValueType>(
                        solutionVectorNonblocked, solutionVectorBlock1, solverTol * 10));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(boundary_operator_agrees_with_2x2_blocked_boundary_operator_for_convergence_testing_in_dual,
                              ValueType, result_types)
{
    typedef ValueType RT;
    typedef typename ScalarTraits<ValueType>::RealType RealType;
    typedef RealType BFT;

    typedef Bempp::DefaultIterativeSolver<BFT, RT> IterSolver;
    const RealType solverTol = 1e-5;

    LaplaceDirichletFixture<BFT, RT> fixture;

    arma::Col<RT> solutionVectorNonblocked;

    // Solve using nonblocked operator
    {
        IterSolver solver(
            fixture.lhsOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);       
        solver.initializeSolver(defaultGmresParameterList(solverTol));
        Solution<BFT, RT> solution = solver.solve(fixture.rhs);
        solutionVectorNonblocked = solution.gridFunction().coefficients();
    }

    // Solve using 2x2 ([A, 0 * A; 0, A]) blocked operator
    {
        BlockedOperatorStructure<BFT, RT> structure;
        structure.setBlock(0, 0, fixture.lhsOp);
        structure.setBlock(0, 1, 0. * fixture.lhsOp);
        structure.setBlock(1, 1, fixture.lhsOp);
        
        BlockedBoundaryOperator<BFT, RT> lhsBlockedOp(structure);
        std::vector<GridFunction<BFT, RT> > blockedRhs(2);
        blockedRhs[0] = fixture.rhs;
        blockedRhs[1] = 2. * fixture.rhs;

        IterSolver solver(
            lhsBlockedOp, IterSolver::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);       
        solver.initializeSolver(defaultGmresParameterList(solverTol));
        BlockedSolution<BFT, RT> solution = solver.solve(blockedRhs);
        arma::Col<RT> solutionVectorBlock0 = solution.gridFunction(0).coefficients();
        arma::Col<RT> solutionVectorBlock1 = solution.gridFunction(1).coefficients() / 2.;

        BOOST_CHECK(check_arrays_are_close<ValueType>(
                        solutionVectorNonblocked, solutionVectorBlock0, solverTol * 10));
        BOOST_CHECK(check_arrays_are_close<ValueType>(
                        solutionVectorNonblocked, solutionVectorBlock1, solverTol * 10));
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
