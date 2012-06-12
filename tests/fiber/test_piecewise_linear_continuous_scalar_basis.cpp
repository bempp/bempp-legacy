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

#include "fiber/piecewise_linear_continuous_scalar_basis.hpp"
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

BOOST_AUTO_TEST_SUITE(PiecewiseLinearContinuousScalarBasis_triangle)

BOOST_AUTO_TEST_CASE_TEMPLATE(size_is_3, ValueType, basis_function_types)
{
    Fiber::PiecewiseLinearContinuousScalarBasis<3, ValueType> basis;
    BOOST_CHECK_EQUAL(basis.size(), 3);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(order_is_1, ValueType, basis_function_types)
{
    Fiber::PiecewiseLinearContinuousScalarBasis<3, ValueType> basis;
    BOOST_CHECK_EQUAL(basis.order(), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, Fiber::ALL_DOFS, data);

    arma::Cube<ValueType> expected(1, // component count
                                   vertexCount,
                                   pointCount);
    expected.fill(0.);
    expected(0, 0, 0) = 1.;
    expected(0, 1, 1) = 1.;
    expected(0, 2, 2) = 1.;

    BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, 1e-10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_second_dof,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = 3;
    const int dofIndex = 1;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, dofIndex, data);

    arma::Cube<ValueType> expected(1, // component count
                                   1, // number of DOFs
                                   pointCount);
    expected.fill(0.);
    expected(0, 0, 1) = 1.;

    BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, 1e-10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_derivatives_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    Fiber::_4dArray<ValueType> expected(1, // component count
                                       elementDim,
                                       vertexCount,
                                       pointCount);
    std::fill(expected.begin(), expected.end(), 0.);
    expected(0, 0, 0, 0) = -1.;
    expected(0, 1, 0, 0) = -1.;
    expected(0, 0, 1, 0) = 1.;
    expected(0, 1, 2, 0) = 1.;
    expected(0, 0, 0, 1) = -1.;
    expected(0, 1, 0, 1) = -1.;
    expected(0, 0, 1, 1) = 1.;
    expected(0, 1, 2, 1) = 1.;
    expected(0, 0, 0, 2) = -1.;
    expected(0, 1, 0, 2) = -1.;
    expected(0, 0, 1, 2) = 1.;
    expected(0, 1, 2, 2) = 1.;

    BOOST_CHECK(check_arrays_are_close<ValueType>(data.derivatives, expected, 1e-10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_and_derivatives_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    const int vertexCount = 3;
    const int elementDim = 2;
    const int pointCount = 3;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES | Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    {
        arma::Cube<ValueType> expected(1, // component count
                                       vertexCount,
                                       pointCount);
        expected.fill(0.);
        expected(0, 0, 0) = 1.;
        expected(0, 1, 1) = 1.;
        expected(0, 2, 2) = 1.;

        BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, 1e-10));
    }
    {
        Fiber::_4dArray<ValueType> expected(1, // component count
                                           elementDim,
                                           vertexCount,
                                           pointCount);
        std::fill(expected.begin(), expected.end(), 0.);
        expected(0, 0, 0, 0) = -1.;
        expected(0, 1, 0, 0) = -1.;
        expected(0, 0, 1, 0) = 1.;
        expected(0, 1, 2, 0) = 1.;
        expected(0, 0, 0, 1) = -1.;
        expected(0, 1, 0, 1) = -1.;
        expected(0, 0, 1, 1) = 1.;
        expected(0, 1, 2, 1) = 1.;
        expected(0, 0, 0, 2) = -1.;
        expected(0, 1, 0, 2) = -1.;
        expected(0, 0, 1, 2) = 1.;
        expected(0, 1, 2, 2) = 1.;
        BOOST_CHECK(check_arrays_are_close<ValueType>(data.derivatives, expected, 1e-10));
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(PiecewiseLinearContinuousScalarBasis_square)

BOOST_AUTO_TEST_CASE_TEMPLATE(size_is_4, ValueType, basis_function_types)
{
    Fiber::PiecewiseLinearContinuousScalarBasis<4, ValueType> basis;
    BOOST_CHECK_EQUAL(basis.size(), 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(order_is_1, ValueType, basis_function_types)
{
    Fiber::PiecewiseLinearContinuousScalarBasis<4, ValueType> basis;
    BOOST_CHECK_EQUAL(basis.order(), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    const int vertexCount = 4;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    points(0, 3) = 1.;
    points(1, 3) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, Fiber::ALL_DOFS, data);

    arma::Cube<ValueType> expected(1, // component count
                                   vertexCount,
                                   pointCount);
    expected.fill(0.);
    expected(0, 0, 0) = 1.;
    expected(0, 1, 1) = 1.;
    expected(0, 2, 2) = 1.;
    expected(0, 3, 3) = 1.;

    BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, 1e-10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_second_dof,
                              ValueType, basis_function_types)
{
    const int vertexCount = 4;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    const int dofIndex = 1;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    points(0, 3) = 1.;
    points(1, 3) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, dofIndex, data);

    arma::Cube<ValueType> expected(1, // component count
                                   1, // number of DOFs
                                   pointCount);
    expected.fill(0.);
    expected(0, 0, 1) = 1.;

    BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, 1e-10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_derivatives_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    const int vertexCount = 4;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    points(0, 3) = 1.;
    points(1, 3) = 1.;

    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    Fiber::_4dArray<ValueType> expected(1, // component count
                                       elementDim,
                                       vertexCount,
                                       pointCount);
    std::fill(expected.begin(), expected.end(), 0.);
    for (int point = 0; point < pointCount; ++point) {
        expected(0, 0, 0, point) = -(1. - points(1, point));
        expected(0, 1, 0, point) = -(1. - points(0, point));
        expected(0, 0, 1, point) = 1. - points(1, point);
        expected(0, 1, 1, point) = -points(0, point);
        expected(0, 0, 2, point) = -points(1, point);
        expected(0, 1, 2, point) = 1. - points(0, point);
        expected(0, 0, 3, point) = points(1, point);
        expected(0, 1, 3, point) = points(0, point);
    }

    BOOST_CHECK(check_arrays_are_close<ValueType>(data.derivatives, expected, 1e-10));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_and_derivatives_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    const int vertexCount = 4;
    const int elementDim = 2;
    const int pointCount = vertexCount;
    typedef Fiber::PiecewiseLinearContinuousScalarBasis<vertexCount, ValueType> Basis;
    Basis basis;
    arma::Mat<typename Basis::CoordinateType> points(elementDim, pointCount);
    points.fill(0.);
    points(0, 1) = 1.;
    points(1, 2) = 1.;
    points(0, 3) = 1.;
    points(1, 3) = 1.;
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES | Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    {
        arma::Cube<ValueType> expected(1, // component count
                                       vertexCount,
                                       pointCount);
        expected.fill(0.);
        expected(0, 0, 0) = 1.;
        expected(0, 1, 1) = 1.;
        expected(0, 2, 2) = 1.;
        expected(0, 3, 3) = 1.;

        BOOST_CHECK(check_arrays_are_close<ValueType>(data.values, expected, 1e-10));
    }
    {
        Fiber::_4dArray<ValueType> expected(1, // component count
                                           elementDim,
                                           vertexCount,
                                           pointCount);
        std::fill(expected.begin(), expected.end(), 0.);
        for (int point = 0; point < pointCount; ++point) {
            expected(0, 0, 0, point) = -(1. - points(1, point));
            expected(0, 1, 0, point) = -(1. - points(0, point));
            expected(0, 0, 1, point) = 1. - points(1, point);
            expected(0, 1, 1, point) = -points(0, point);
            expected(0, 0, 2, point) = -points(1, point);
            expected(0, 1, 2, point) = 1. - points(0, point);
            expected(0, 0, 3, point) = points(1, point);
            expected(0, 1, 3, point) = points(0, point);
        }

        BOOST_CHECK(check_arrays_are_close<ValueType>(data.derivatives, expected, 1e-10));
    }
}

BOOST_AUTO_TEST_SUITE_END()
