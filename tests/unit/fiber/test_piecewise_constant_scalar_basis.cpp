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

#include "fiber/piecewise_constant_scalar_basis.hpp"
#include "../type_template.hpp"

#include <algorithm>
#include "common/eigen_support.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>

// Tests

using namespace Bempp;

BOOST_AUTO_TEST_SUITE(PiecewiseConstantScalarBasis)

BOOST_AUTO_TEST_CASE_TEMPLATE(size_is_1, ValueType, basis_function_types)
{
    Fiber::PiecewiseConstantScalarBasis<ValueType> basis;
    BOOST_CHECK_EQUAL(basis.size(), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(order_is_0, ValueType, basis_function_types)
{
    Fiber::PiecewiseConstantScalarBasis<ValueType> basis;
    BOOST_CHECK_EQUAL(basis.order(), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    typedef Fiber::PiecewiseConstantScalarBasis<ValueType> Basis;
    Basis basis;
    Matrix<typename Basis::CoordinateType> points(3, 4);
    srand(1);
    points.setRandom();
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, Fiber::ALL_DOFS, data);

    Fiber::_3dArray<ValueType> expected(1, 1, points.cols());
    std::fill(expected.begin(), expected.end(), 1.);

    BOOST_CHECK_EQUAL(data.values.extent(0), 1u); // 1 component
    BOOST_CHECK_EQUAL(data.values.extent(1), 1u); // 1 basis function
    BOOST_CHECK_EQUAL(data.values.extent(2), points.cols());
    BOOST_CHECK_EQUAL_COLLECTIONS(data.values.begin(),
                                  data.values.end(),
                                  expected.begin(),
                                  expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_works_for_one_dof,
                              ValueType, basis_function_types)
{
    typedef Fiber::PiecewiseConstantScalarBasis<ValueType> Basis;
    Basis basis;
    Matrix<typename Basis::CoordinateType> points(3, 4);
    srand(1);
    points.setRandom();
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES, points, 0 /* dof number */, data);

    Fiber::_3dArray<ValueType> expected(1, 1, points.cols());
    std::fill(expected.begin(), expected.end(), 1.);

    BOOST_CHECK_EQUAL(data.values.extent(0), 1u); // 1 component
    BOOST_CHECK_EQUAL(data.values.extent(1), 1u); // 1 basis function
    BOOST_CHECK_EQUAL(data.values.extent(2), points.cols());
    BOOST_CHECK_EQUAL_COLLECTIONS(data.values.begin(),
                                  data.values.end(),
                                  expected.begin(),
                                  expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_derivatives_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    typedef Fiber::PiecewiseConstantScalarBasis<ValueType> Basis;
    Basis basis;
    Matrix<typename Basis::CoordinateType> points(3, 4);
    srand(1);
    points.setRandom();
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::DERIVATIVES, points, Fiber::ALL_DOFS, data);

    Fiber::_4dArray<ValueType> expected(1,  // 1 component
                                        points.rows(),
                                        1,  // 1 basis functions
                                        points.cols());
    std::fill(expected.begin(), expected.end(), 0.);
    BOOST_CHECK_EQUAL(data.derivatives.extent(0), expected.extent(0));
    BOOST_CHECK_EQUAL(data.derivatives.extent(1), expected.extent(1));
    BOOST_CHECK_EQUAL(data.derivatives.extent(2), expected.extent(2));
    BOOST_CHECK_EQUAL(data.derivatives.extent(3), expected.extent(3));
    BOOST_CHECK_EQUAL_COLLECTIONS(data.derivatives.begin(),
                                  data.derivatives.end(),
                                  expected.begin(),
                                  expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_derivatives_works_for_one_dof,
                              ValueType, basis_function_types)
{
    typedef Fiber::PiecewiseConstantScalarBasis<ValueType> Basis;
    Basis basis;
    Matrix<typename Basis::CoordinateType> points(3, 4);
    srand(1);
    points.setRandom();
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::DERIVATIVES, points, 0 /* dof number */, data);

    Fiber::_4dArray<ValueType> expected(1,  // 1 component
                                        points.rows(),
                                        1,  // 1 basis functions
                                        points.cols());
    std::fill(expected.begin(), expected.end(), 0.);
    BOOST_CHECK_EQUAL(data.derivatives.extent(0), expected.extent(0));
    BOOST_CHECK_EQUAL(data.derivatives.extent(1), expected.extent(1));
    BOOST_CHECK_EQUAL(data.derivatives.extent(2), expected.extent(2));
    BOOST_CHECK_EQUAL(data.derivatives.extent(3), expected.extent(3));
    BOOST_CHECK_EQUAL_COLLECTIONS(data.derivatives.begin(),
                                  data.derivatives.end(),
                                  expected.begin(),
                                  expected.end());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluate_values_and_derivatives_works_for_all_dofs,
                              ValueType, basis_function_types)
{
    typedef Fiber::PiecewiseConstantScalarBasis<ValueType> Basis;
    Basis basis;
    Matrix<typename Basis::CoordinateType> points(3, 4);
    srand(1);
    points.setRandom();
    Fiber::BasisData<ValueType> data;
    basis.evaluate(Fiber::VALUES | Fiber::DERIVATIVES, points, 0 /* dof number */, data);

    {
        Fiber::_3dArray<ValueType> expected(1, 1, points.cols());
        std::fill(expected.begin(), expected.end(), 1.);
        BOOST_CHECK_EQUAL(data.values.extent(0), 1); // 1 component
        BOOST_CHECK_EQUAL(data.values.extent(1), 1); // 1 basis function
        BOOST_CHECK_EQUAL(data.values.extent(2), points.cols());
        BOOST_CHECK_EQUAL_COLLECTIONS(data.values.begin(),
                                      data.values.end(),
                                      expected.begin(),
                                      expected.end());
    }

    {
        Fiber::_4dArray<ValueType> expected(1,  // 1 component
                                            points.rows(),
                                            1,  // 1 basis functions
                                            points.cols());
        std::fill(expected.begin(), expected.end(), 0.);
        BOOST_CHECK_EQUAL(data.derivatives.extent(0), expected.extent(0));
        BOOST_CHECK_EQUAL(data.derivatives.extent(1), expected.extent(1));
        BOOST_CHECK_EQUAL(data.derivatives.extent(2), expected.extent(2));
        BOOST_CHECK_EQUAL(data.derivatives.extent(3), expected.extent(3));
        BOOST_CHECK_EQUAL_COLLECTIONS(data.derivatives.begin(),
                                      data.derivatives.end(),
                                      expected.begin(),
                                      expected.end());
    }
}

BOOST_AUTO_TEST_SUITE_END()
