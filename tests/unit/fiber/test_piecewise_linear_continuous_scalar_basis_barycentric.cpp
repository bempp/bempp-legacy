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

#include "fiber/linear_scalar_shapeset_barycentric.hpp"
#include "fiber/scalar_traits.hpp"
#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

BOOST_AUTO_TEST_SUITE(LinearScalarShapesetBarycentric_triangle)

BOOST_AUTO_TEST_CASE_TEMPLATE(size_is_3, ValueType, basis_function_types)
{
    typedef Fiber::LinearScalarShapesetBarycentric<ValueType> Shapeset;
    Shapeset basis(Shapeset::TYPE1);
    BOOST_CHECK_EQUAL(basis.size(), 3);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(order_is_1, ValueType, basis_function_types)
{
    typedef Fiber::LinearScalarShapesetBarycentric<ValueType> Shapeset;
    Shapeset basis(Shapeset::TYPE1);
    BOOST_CHECK_EQUAL(basis.order(), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_all_dofs_type1,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::LinearScalarShapesetBarycentric<ValueType> Shapeset;
    Shapeset basis(Shapeset::TYPE1);
    typedef typename Shapeset::CoordinateType CoordinateType;
    arma::Mat<CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, Fiber::ALL_DOFS, data);

    Fiber::_3dArray<ValueType> expected(1, // component count
                                        vertexCount,
                                        pointCount);
    std::fill(expected.begin(), expected.end(), 0.);
    expected(0, 0, 0) = 1.;
    expected(0, 0, 1) = 1./3;
    expected(0, 0, 2) = 1./2;

    expected(0, 1, 0) = 0;
    expected(0, 1, 1) = 1./3;
    expected(0, 1, 2) = 0;

    expected(0, 2, 0) = 0;
    expected(0, 2, 1) = 1./3;
    expected(0, 2, 2) = 1./2;

    const CoordinateType tol =
        10 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, tol));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_all_dofs_type2,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::LinearScalarShapesetBarycentric<ValueType> Shapeset;
    Shapeset basis(Shapeset::TYPE2);
    typedef typename Shapeset::CoordinateType CoordinateType;
    arma::Mat<CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, Fiber::ALL_DOFS, data);

    Fiber::_3dArray<ValueType> expected(1, // component count
                                        vertexCount,
                                        pointCount);
    std::fill(expected.begin(), expected.end(), 0.);
    expected(0, 0, 0) = 1.;
    expected(0, 0, 1) = 1./2;
    expected(0, 0, 2) = 1./3;

    expected(0, 1, 0) = 0;
    expected(0, 1, 1) = 1./2;
    expected(0, 1, 2) = 1./3;

    expected(0, 2, 0) = 0;
    expected(0, 2, 1) = 0;
    expected(0, 2, 2) = 1./3;

    const CoordinateType tol =
        10 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, tol));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_derivatives_works_for_all_dofs_type1,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::LinearScalarShapesetBarycentric<ValueType> Shapeset;
    typedef typename Shapeset::CoordinateType CoordinateType;
    Shapeset basis(Shapeset::TYPE1);
    arma::Mat<CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    Fiber::_4dArray<ValueType> expected(1, // component count
                                        2, //
                                        vertexCount,
                                        pointCount);
    std::fill(expected.begin(), expected.end(), 0.);

    expected(0, 0, 0, 0) = -2./3.;
    expected(0, 1, 0, 0) = -1./2.;
    expected(0, 0, 0, 1) = -2./3.;
    expected(0, 1, 0, 1) = -1./2.;
    expected(0, 0, 0, 2) = -2./3.;
    expected(0, 1, 0, 2) = -1./2.;

    expected(0, 0, 1, 0) = 1./3.;
    expected(0, 1, 1, 0) = 0.;
    expected(0, 0, 1, 1) = 1./3.;
    expected(0, 1, 1, 1) = 0.;
    expected(0, 0, 1, 2) = 1./3.;
    expected(0, 1, 1, 2) = 0.;

    expected(0, 0, 2, 0) = 1./3.;
    expected(0, 1, 2, 0) = 1./2;
    expected(0, 0, 2, 1) = 1./3.;
    expected(0, 1, 2, 1) = 1./2;
    expected(0, 0, 2, 2) = 1./3.;
    expected(0, 1, 2, 2) = 1./2;

    const CoordinateType tol =
        10 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(data.derivatives, expected, 
						  tol));
}


BOOST_AUTO_TEST_SUITE_END()

