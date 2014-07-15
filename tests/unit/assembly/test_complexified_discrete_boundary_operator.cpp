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

template <typename RealType>
struct ComplexifiedDiscreteBoundaryOperatorFixture
{
    ComplexifiedDiscreteBoundaryOperatorFixture()
    {
        arma::Mat<RealType> mat = generateRandomMatrix<RealType>(4, 5);
        op.reset(new DiscreteDenseBoundaryOperator<RealType>(mat));
        complexifiedOp = complexify(op);
    }

    shared_ptr<const DiscreteBoundaryOperator<RealType> > op;
    shared_ptr<const DiscreteBoundaryOperator<std::complex<RealType> > > complexifiedOp;
};

} // namespace

BOOST_AUTO_TEST_SUITE(ComplexifiedDiscreteBoundaryOperator)

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_alpha_equal_to_2_plus_3j_and_beta_equal_to_0_and_y_initialized_to_nans,
                              RealType, real_result_types)
{
    std::srand(1);

    typedef std::complex<RealType> ComplexType;
    ComplexType complexNan(std::numeric_limits<RealType>::quiet_NaN(),
                           std::numeric_limits<RealType>::quiet_NaN());

    ComplexifiedDiscreteBoundaryOperatorFixture<RealType> fixture;
    arma::Mat<RealType> mat = fixture.op->asMatrix();
    arma::Mat<ComplexType> complexMat(mat.n_rows, mat.n_cols);
    complexMat.fill(0.);
    complexMat.set_real(mat);

    shared_ptr<const DiscreteBoundaryOperator<ComplexType> > dop = fixture.complexifiedOp;

    ComplexType alpha(2., 3.);
    ComplexType beta(0.);

    arma::Col<ComplexType> x = generateRandomVector<ComplexType>(dop->columnCount());
    arma::Col<ComplexType> y(dop->rowCount());
    y.fill(complexNan);

    arma::Col<ComplexType> expected = alpha * mat * x;
    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(y.is_finite());
    BOOST_CHECK(check_arrays_are_close<ComplexType>(
                    y, expected,  10. * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(builtin_apply_works_correctly_for_no_transpose_andalpha_equal_to_2_plus_3j_and_beta_equal_to_4_minus_5j,
                              RealType, real_result_types)
{
    std::srand(1);

    typedef std::complex<RealType> ComplexType;

    ComplexifiedDiscreteBoundaryOperatorFixture<RealType> fixture;
    arma::Mat<RealType> mat = fixture.op->asMatrix();
    arma::Mat<ComplexType> complexMat(mat.n_rows, mat.n_cols);
    complexMat.fill(0.);
    complexMat.set_real(mat);

    shared_ptr<const DiscreteBoundaryOperator<ComplexType> > dop = fixture.complexifiedOp;

    ComplexType alpha(2., 3.);
    ComplexType beta(4., -5.);

    arma::Col<ComplexType> x = generateRandomVector<ComplexType>(dop->columnCount());
    arma::Col<ComplexType> y = generateRandomVector<ComplexType>(dop->rowCount());

    arma::Col<ComplexType> expected = alpha * mat * x + beta * y;

    dop->apply(NO_TRANSPOSE, x, y, alpha, beta);
    
    BOOST_CHECK(check_arrays_are_close<ComplexType>(
                    y, expected, 10. * std::numeric_limits<RealType>::epsilon()));
}

BOOST_AUTO_TEST_SUITE_END()
