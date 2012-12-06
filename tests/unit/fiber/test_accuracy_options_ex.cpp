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

#include "fiber/accuracy_options.hpp"
#include "fiber/quadrature_options.hpp"

#include <boost/test/unit_test.hpp>

// Tests

using namespace Fiber;

BOOST_AUTO_TEST_SUITE(AccuracyOptionsEx)

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_uniform_relative_order)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    opts.setDoubleRegular(order1);
    int order = opts.doubleRegular(2.).quadratureOrder(defaultOrder);

    BOOST_CHECK_EQUAL(order, defaultOrder + order1);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_uniform_absolute_order)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    opts.setDoubleRegular(order1, false /* absolute */);
    int order = opts.doubleRegular(2.).quadratureOrder(defaultOrder);

    BOOST_CHECK_EQUAL(order, order1);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_two_relative_orders)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    const double maxNormalizedDistance1 = 4.;
    const int order2 = 2;
    opts.setDoubleRegular(maxNormalizedDistance1, order1, order2);

    int orderNear = opts.doubleRegular(maxNormalizedDistance1 - 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderNear, defaultOrder + order1);

    int orderFar = opts.doubleRegular(maxNormalizedDistance1 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderFar, defaultOrder + order2);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_three_absolute_orders)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    const double maxNormalizedDistance1 = 5.;
    const int order2 = 2;
    const double maxNormalizedDistance2 = 10.;
    const int order3 = 1;
    opts.setDoubleRegular(maxNormalizedDistance1, order1,
                          maxNormalizedDistance2, order2,
                          order3, false /* absolute */);

    int orderNear = opts.doubleRegular(maxNormalizedDistance1 - 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderNear, order1);

    int orderMiddle = opts.doubleRegular(maxNormalizedDistance1 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle, order2);

    int orderFar = opts.doubleRegular(maxNormalizedDistance2 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderFar, order3);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_four_relative_orders)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    const double maxNormalizedDistance1 = 5.;
    const int order2 = 2;
    const double maxNormalizedDistance2 = 10.;
    const int order3 = 1;
    const double maxNormalizedDistance3 = 15.;
    const int order4 = 0;
    opts.setDoubleRegular(maxNormalizedDistance1, order1,
                          maxNormalizedDistance2, order2,
                          maxNormalizedDistance3, order3,
                          order4);

    int orderNear = opts.doubleRegular(maxNormalizedDistance1 - 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderNear, defaultOrder + order1);

    int orderMiddle1 = opts.doubleRegular(maxNormalizedDistance1 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle1, defaultOrder + order2);

    int orderMiddle2 = opts.doubleRegular(maxNormalizedDistance2 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle2, defaultOrder + order3);

    int orderFar = opts.doubleRegular(maxNormalizedDistance3 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderFar, defaultOrder + order4);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_five_absolute_orders)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    const double maxNormalizedDistance1 = 5.;
    const int order2 = 2;
    const double maxNormalizedDistance2 = 10.;
    const int order3 = 1;
    const double maxNormalizedDistance3 = 15.;
    const int order4 = 0;
    const double maxNormalizedDistance4 = 20.;
    const int order5 = -1;
    opts.setDoubleRegular(maxNormalizedDistance1, order1,
                          maxNormalizedDistance2, order2,
                          maxNormalizedDistance3, order3,
                          maxNormalizedDistance4, order4,
                          order5, false /* absolute */);

    int orderNear = opts.doubleRegular(maxNormalizedDistance1 - 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderNear, order1);

    int orderMiddle1 = opts.doubleRegular(maxNormalizedDistance1 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle1, order2);

    int orderMiddle2 = opts.doubleRegular(maxNormalizedDistance2 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle2, order3);

    int orderMiddle3 = opts.doubleRegular(maxNormalizedDistance3 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle3, order4);

    int orderFar = opts.doubleRegular(maxNormalizedDistance4 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderFar, order5);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_three_absolute_orders_given_as_vectors)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    const double maxNormalizedDistance1 = 5.;
    const int order2 = 2;
    const double maxNormalizedDistance2 = 10.;
    const int order3 = 1;
    std::vector<double> maxNormalizedDistances;
    maxNormalizedDistances.push_back(maxNormalizedDistance1);
    maxNormalizedDistances.push_back(maxNormalizedDistance2);
    std::vector<int> orders;
    orders.push_back(order1);
    orders.push_back(order2);
    orders.push_back(order3);

    opts.setDoubleRegular(maxNormalizedDistances, orders, false /* absolute */);

    int orderNear = opts.doubleRegular(maxNormalizedDistance1 - 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderNear, order1);

    int orderMiddle = opts.doubleRegular(maxNormalizedDistance1 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle, order2);

    int orderFar = opts.doubleRegular(maxNormalizedDistance2 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderFar, order3);
}

BOOST_AUTO_TEST_CASE(quadratureOrder_agrees_with_setDoubleRegular_for_three_relative_orders_given_as_vectors)
{
    Fiber::AccuracyOptionsEx opts;
    const int defaultOrder = 3;
    const int order1 = 5;
    const double maxNormalizedDistance1 = 5.;
    const int order2 = 2;
    const double maxNormalizedDistance2 = 10.;
    const int order3 = 1;
    std::vector<double> maxNormalizedDistances;
    maxNormalizedDistances.push_back(maxNormalizedDistance1);
    maxNormalizedDistances.push_back(maxNormalizedDistance2);
    std::vector<int> orders;
    orders.push_back(order1);
    orders.push_back(order2);
    orders.push_back(order3);

    opts.setDoubleRegular(maxNormalizedDistances, orders);

    int orderNear = opts.doubleRegular(maxNormalizedDistance1 - 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderNear, defaultOrder + order1);

    int orderMiddle = opts.doubleRegular(maxNormalizedDistance1 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderMiddle, defaultOrder + order2);

    int orderFar = opts.doubleRegular(maxNormalizedDistance2 + 0.1)
            .quadratureOrder(defaultOrder);
    BOOST_CHECK_EQUAL(orderFar, defaultOrder + order3);
}

BOOST_AUTO_TEST_SUITE_END()
