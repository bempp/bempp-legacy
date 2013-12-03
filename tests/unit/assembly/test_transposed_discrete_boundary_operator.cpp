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

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"
#include "../random_arrays.hpp"

#include "create_regular_grid.hpp"

#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/discrete_dense_boundary_operator.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

using namespace Bempp;

namespace 
{

template <typename RT>
struct TransposedDiscreteBoundaryOperatorFixture
{
    TransposedDiscreteBoundaryOperatorFixture()
    {
        arma::Mat<RT> mat = generateRandomMatrix<RT>(4, 5);
        op.reset(new DiscreteDenseBoundaryOperator<RT>(mat));
        unmodifiedOp = transpose(NO_TRANSPOSE, op);
        transposedOp = transpose(op);
        conjugatedOp = conjugate(op);
        conjugateTransposedOp = conjugateTranspose(op);
    }

    shared_ptr<const DiscreteBoundaryOperator<RT> >
    op, unmodifiedOp, transposedOp, conjugatedOp, conjugateTransposedOp;
};

} // namespace

BOOST_AUTO_TEST_SUITE(TransposedDiscreteBoundaryOperator)

// UNMODIFIED //////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_no_transpose_and_alpha_equal_to_2_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.unmodifiedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_no_transpose_and_alpha_equal_to_2_plus_3j_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.unmodifiedOp;

    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_no_transpose_and_alpha_equal_to_2_and_beta_equal_to_3,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.unmodifiedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_no_transpose_andalpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.unmodifiedOp;

    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    std::complex<double> tempb(4.,-5.);
    RT beta = static_cast<RT>(tempb);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected, 
                                           10. * std::numeric_limits<CT>::epsilon()));
}

// TRANSPOSED //////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_transpose_and_alpha_equal_to_2_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().st();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.transposedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_transpose_and_alpha_equal_to_2_plus_3j_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().st();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.transposedOp;

    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_transpose_and_alpha_equal_to_2_and_beta_equal_to_3,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().st();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.transposedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_transpose_andalpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().st();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.transposedOp;

    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    std::complex<double> tempb(4.,-5.);
    RT beta = static_cast<RT>(tempb);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

// CONJUGATED //////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_and_alpha_equal_to_2_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = arma::conj(fixture.op->asMatrix());
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugatedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_and_alpha_equal_to_2_plus_3j_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = arma::conj(fixture.op->asMatrix());
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugatedOp;


    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_and_alpha_equal_to_2_and_beta_equal_to_3,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = arma::conj(fixture.op->asMatrix());
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugatedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_andalpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = arma::conj(fixture.op->asMatrix());
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugatedOp;

    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    std::complex<double> tempb(4.,-5.);
    RT beta = static_cast<RT>(tempb);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

// CONJUGATE_TRANSPOSED ////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_transpose_and_alpha_equal_to_2_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().t();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugateTransposedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_transpose_and_alpha_equal_to_2_plus_3j_and_beta_equal_to_0_and_y_initialized_to_nans,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().t();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugateTransposedOp;

    RT alpha = static_cast<RT>(2., 3.);
    RT beta = static_cast<RT>(0.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y(dop->rowCount());
    y.fill(std::numeric_limits<CT>::quiet_NaN());

    arma::Col<RT> expected = alpha * mat * x;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_transpose_and_alpha_equal_to_2_and_beta_equal_to_3,
                              ResultType, result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().t();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugateTransposedOp;

    RT alpha = static_cast<RT>(2.);
    RT beta = static_cast<RT>(3.);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_conjugate_transpose_andalpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j,
                              ResultType, complex_result_types)
{
    std::srand(1);

    typedef ResultType RT;
    typedef typename Fiber::ScalarTraits<RT>::RealType BFT;
    typedef typename Fiber::ScalarTraits<RT>::RealType CT;

    TransposedDiscreteBoundaryOperatorFixture<RT> fixture;
    arma::Mat<RT> mat = fixture.op->asMatrix().t();
    shared_ptr<const DiscreteBoundaryOperator<RT> > dop = fixture.conjugateTransposedOp;

    std::complex<double> tempa(2.,3.);
    RT alpha = static_cast<RT>(tempa);
    std::complex<double> tempb(4.,-5.);
    RT beta = static_cast<RT>(tempb);

    arma::Col<RT> x = generateRandomVector<RT>(dop->columnCount());
    arma::Col<RT> y = generateRandomVector<RT>(dop->rowCount());

    arma::Col<RT> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);

    BOOST_CHECK(check_arrays_are_close<RT>(y, expected,
                                           10. * std::numeric_limits<CT>::epsilon()));
}

BOOST_AUTO_TEST_SUITE_END()
