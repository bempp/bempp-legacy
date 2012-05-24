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

#include "assembly/discrete_dense_linear_operator.hpp"
#include "linalg/belos_solver_wrapper.hpp"

#include <armadillo>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/version.hpp>

#include <Thyra_DefaultSpmdVector.hpp>

// Real ValueType
template <typename ValueType>
typename boost::disable_if<
boost::is_complex<ValueType>,
arma::Mat<ValueType>
>::type
generateRandomMatrix(int rowCount, int colCount)
{
    return arma::randu<arma::Mat<ValueType> >(rowCount, colCount);
}

template <typename ValueType>
typename boost::enable_if<
boost::is_complex<ValueType>,
arma::Mat<ValueType>
>::type
generateRandomMatrix(int rowCount, int colCount)
{
    return arma::randu<arma::Mat<ValueType> >(rowCount, colCount) +
            ValueType(0, 1) *
            arma::randu<arma::Mat<ValueType> >(rowCount, colCount);
}

// Tests

BOOST_AUTO_TEST_SUITE(BelosSolverWrapper)

// THIS TEST FAILS FOR float AND std::complex<float> DUE TO A BUG IN TRILINOS.
// I DON'T KNOW YET HOW TO REGISTER EXPECTED FAILURES FOR TEMPLATED TESTS.

BOOST_AUTO_TEST_CASE_TEMPLATE(solve_works_for_single_rhs, ValueType, result_types)
{
    std::srand(1);

    const int size = 50;
    arma::Mat<ValueType> mat = generateRandomMatrix<ValueType>(size, size);
    Bempp::DiscreteDenseLinearOperator<ValueType> op(mat);

    arma::Col<ValueType> rhs = generateRandomMatrix<ValueType>(size, 1);
    arma::Col<ValueType> sol(size);
    sol.fill(static_cast<ValueType>(0.));

    typedef Thyra::DefaultSpmdVector<ValueType> DenseVector;
    Teuchos::ArrayRCP<ValueType> rhsArray =
            Teuchos::arcp(rhs.memptr(), 0 /* lowerOffset */,
                          size, false /* doesn't own memory */);
    DenseVector rhsVector(Thyra::defaultSpmdVectorSpace<ValueType>(size),
                    rhsArray, 1 /* stride */);
    Teuchos::ArrayRCP<ValueType> solArray =
            Teuchos::arcp(sol.memptr(), 0 /* lowerOffset */,
                          size, false /* doesn't own memory */);
    DenseVector solVector(Thyra::defaultSpmdVectorSpace<ValueType>(size),
                    solArray, 1 /* stride */);

    typedef Bempp::BelosSolverWrapper<ValueType> Solver;
    Solver solver(Teuchos::rcpFromRef<const Thyra::LinearOpBase<ValueType> >(op));
    solver.initializeSolver(Bempp::defaultGmresParameterList(1e-10));

    Thyra::SolveStatus<typename Solver::MagnitudeType > status =
            solver.solve(Thyra::NOTRANS, rhsVector,
                         Teuchos::ptr<Thyra::MultiVectorBase<ValueType> >(&solVector));
    BOOST_CHECK_EQUAL(status.solveStatus, Thyra::SOLVE_STATUS_CONVERGED);

    // Solve system with armadillo

    arma::Col<ValueType> armaSol = arma::solve(mat, rhs);

    BOOST_CHECK(check_arrays_are_close<ValueType>(sol, armaSol, 1e-6));
}

BOOST_AUTO_TEST_SUITE_END()

#endif
