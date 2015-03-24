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

#include "bempp/common/config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"
#include "../random_arrays.hpp"

#include "create_regular_grid.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/discrete_boundary_operator_composition.hpp"
#include "assembly/discrete_sparse_boundary_operator.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/identity_operator.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/sparse_cholesky.hpp"

#include "bempp/common/config_ahmed.hpp"

#include "grid/grid.hpp"

#include "space/piecewise_constant_scalar_space.hpp"
#include "space/piecewise_linear_discontinuous_scalar_space.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include "common/eigen_support.hpp"
#include "common/boost_make_shared_fwd.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

using namespace Bempp;

// Tests

BOOST_AUTO_TEST_SUITE(SparseCholesky)

BOOST_AUTO_TEST_CASE(sparse_cholesky_works_for_piecewise_linears)
{
    typedef double BFT;
    typedef double RT;

    int nElementsX = 1, nElementsY = 2;
    shared_ptr<Grid> grid = createRegularTriangularGrid(nElementsX, nElementsY);

    shared_ptr<Space<BFT> > space(
        new PiecewiseLinearDiscontinuousScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    BoundaryOperator<BFT, RT> op = identityOperator<BFT, RT>(
        context, space, space, space);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();
    typedef DiscreteSparseBoundaryOperator<RT> SparseOp;
    shared_ptr<const SparseOp> sdop = SparseOp::castToSparse(dop);
    shared_ptr<const Epetra_CrsMatrix> mat = sdop->epetraMatrix();
    shared_ptr<Epetra_CrsMatrix> L = sparseCholesky(*mat);

    shared_ptr<SparseOp> opL = boost::make_shared<SparseOp>(L, NO_SYMMETRY);
    shared_ptr<SparseOp> opLT = boost::make_shared<SparseOp>(L, NO_SYMMETRY,
                                                             TRANSPOSE);
    typedef DiscreteBoundaryOperatorComposition<RT> OpComposition;
    shared_ptr<DiscreteBoundaryOperatorComposition<RT> > opLLT =
            boost::make_shared<OpComposition>(opL, opLT);

    Matrix<double> A = dop->asMatrix();
    Matrix<double> LLT = opLLT->asMatrix();

//    arma::Mat<double> matL = opL->asMatrix();
//    std::cout << "A:\n" << A;
//    std::cout << "L:\n" << matL;
//    std::cout << "arma chol:\n" << arma::chol(A).t();
//    std::cout << "L*L.t():\n" << (matL * matL.t());
//    std::cout << "LLT:\n" << LLT;

    BOOST_CHECK(check_arrays_are_close<double>(
                    A, LLT, 100. * std::numeric_limits<double>::epsilon()));
}

BOOST_AUTO_TEST_CASE(sparse_cholesky_works_for_piecewise_constants)
{
    typedef double BFT;
    typedef double RT;

    int nElementsX = 1, nElementsY = 2;
    shared_ptr<Grid> grid = createRegularTriangularGrid(nElementsX, nElementsY);

    shared_ptr<Space<BFT> > space(
        new PiecewiseConstantScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    BoundaryOperator<BFT, RT> op = identityOperator<BFT, RT>(
        context, space, space, space);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();
    typedef DiscreteSparseBoundaryOperator<RT> SparseOp;
    shared_ptr<const SparseOp> sdop = SparseOp::castToSparse(dop);
    shared_ptr<const Epetra_CrsMatrix> mat = sdop->epetraMatrix();
    shared_ptr<Epetra_CrsMatrix> L = sparseCholesky(*mat);

    shared_ptr<SparseOp> opL = boost::make_shared<SparseOp>(L, NO_SYMMETRY);
    shared_ptr<SparseOp> opLT = boost::make_shared<SparseOp>(L, NO_SYMMETRY,
                                                             TRANSPOSE);
    typedef DiscreteBoundaryOperatorComposition<RT> OpComposition;
    shared_ptr<DiscreteBoundaryOperatorComposition<RT> > opLLT =
            boost::make_shared<OpComposition>(opL, opLT);

    Matrix<double> A = dop->asMatrix();
    Matrix<double> LLT = opLLT->asMatrix();

    BOOST_CHECK(check_arrays_are_close<double>(
                    A, LLT, 100. * std::numeric_limits<double>::epsilon()));
}

BOOST_AUTO_TEST_CASE(sparse_cholesky_throws_for_continuous_functions)
{
    typedef double BFT;
    typedef double RT;

    int nElementsX = 1, nElementsY = 2;
    shared_ptr<Grid> grid = createRegularTriangularGrid(nElementsX, nElementsY);

    shared_ptr<Space<BFT> > space(
        new PiecewiseLinearContinuousScalarSpace<BFT>(grid));

    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    BoundaryOperator<BFT, RT> op = identityOperator<BFT, RT>(
        context, space, space, space);
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = op.weakForm();
    typedef DiscreteSparseBoundaryOperator<RT> SparseOp;
    shared_ptr<const SparseOp> sdop = SparseOp::castToSparse(dop);
    shared_ptr<const Epetra_CrsMatrix> mat = sdop->epetraMatrix();
    BOOST_CHECK_THROW(sparseCholesky(*mat), std::invalid_argument   );
}

BOOST_AUTO_TEST_SUITE_END()

#endif // WITH_TRILINOS
