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

#include "fiber/laplace_3d_single_layer_potential_kernel.hpp"
#include "fiber/geometrical_data.hpp"
#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"

#include <algorithm>
#include <armadillo>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

BOOST_AUTO_TEST_SUITE(Laplace3dSingleLayerPotentialKernel)

BOOST_AUTO_TEST_CASE_TEMPLATE(worldDimension_is_3, ValueType, kernel_types)
{
    Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> op;
    size_t val=3;
    BOOST_CHECK_EQUAL(op.worldDimension(), val);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(domainDimension_is_1, ValueType, kernel_types)
{
    Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> op;
    size_t val=1;
    BOOST_CHECK_EQUAL(op.domainDimension(), val);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(codomainDimension_is_1, ValueType, kernel_types)
{
    Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> op;
    size_t val=1;
    BOOST_CHECK_EQUAL(op.domainDimension(), val);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(addGeometricalDependencies_works,
                              ValueType, kernel_types)
{
    Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> op;
    size_t testGeomDeps = 1024, trialGeomDeps = 16; // random initial values
    op.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    BOOST_CHECK(testGeomDeps & Fiber::GLOBALS);
    BOOST_CHECK(trialGeomDeps & Fiber::GLOBALS);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_works_for_points_on_x_axis,
                              ValueType, kernel_types)
{
    typedef Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> Operator;
    Operator op;
    Fiber::GeometricalData<typename Operator::CoordinateType> testGeomData, trialGeomData;
    const size_t worldDim = 3;
    const size_t testPointCount = 2, trialPointCount = 3;
    // Note: points lying on x axis only
    testGeomData.globals.set_size(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 1) = 1.;

    trialGeomData.globals.set_size(worldDim, trialPointCount);
    trialGeomData.globals.fill(0.);
    trialGeomData.globals(0, 0) = 2.;
    trialGeomData.globals(0, 1) = 3.;
    trialGeomData.globals(0, 2) = 4.;

    Fiber::Array4d<ValueType> result;
    op.evaluateOnGrid(testGeomData, trialGeomData, result);

    Fiber::Array4d<ValueType> expected(1, testPointCount, 1, trialPointCount);
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
        for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint)
            expected(0, testPoint, 0, trialPoint) =
                    (1. / (4. * M_PI)) /
                    std::abs(testGeomData.globals(0, testPoint) -
                             trialGeomData.globals(0, trialPoint));

    BOOST_CHECK(check_arrays_are_close<ValueType>(result, expected, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_works_for_points_on_y_axis,
                              ValueType, kernel_types)
{
    typedef Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> Operator;
    Operator op;
    Fiber::GeometricalData<typename Operator::CoordinateType> testGeomData, trialGeomData;
    const size_t worldDim = 3;
    const size_t testPointCount = 2, trialPointCount = 3;
    // Note: points lying on y axis only
    testGeomData.globals.set_size(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(1, 1) = 1.;

    trialGeomData.globals.set_size(worldDim, trialPointCount);
    trialGeomData.globals.fill(0.);
    trialGeomData.globals(1, 0) = 2.;
    trialGeomData.globals(1, 1) = 3.;
    trialGeomData.globals(1, 2) = 4.;

    Fiber::Array4d<ValueType> result;
    op.evaluateOnGrid(testGeomData, trialGeomData, result);

    Fiber::Array4d<ValueType> expected(1, testPointCount, 1, trialPointCount);
    for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
        for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint)
            expected(0, testPoint, 0, trialPoint) =
                    (1. / (4. * M_PI)) /
                    std::abs(testGeomData.globals(0, testPoint) -
                             trialGeomData.globals(0, trialPoint));

    BOOST_CHECK(check_arrays_are_close<ValueType>(result, expected, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_agrees_with_evaluateAtPointPairs,
                              ValueType, kernel_types)
{
    typedef Fiber::Laplace3dSingleLayerPotentialKernel<ValueType> Operator;
    Operator op;

    typedef Fiber::GeometricalData<typename Operator::CoordinateType> GeomData;

    const size_t worldDim = 3;
    const size_t testPointCount = 2, trialPointCount = 3;

    // Collect data with evaluateOnGrid
    GeomData testGeomDataOnGrid, trialGeomDataOnGrid;

    testGeomDataOnGrid.globals.set_size(worldDim, testPointCount);
    testGeomDataOnGrid.globals.fill(0.);
    testGeomDataOnGrid.globals(0, 1) = 1.;

    trialGeomDataOnGrid.globals.set_size(worldDim, trialPointCount);
    trialGeomDataOnGrid.globals.fill(1.);
    trialGeomDataOnGrid.globals(0, 0) = 2.;
    trialGeomDataOnGrid.globals(0, 1) = 3.;
    trialGeomDataOnGrid.globals(0, 2) = 4.;

    Fiber::Array4d<ValueType> resultOnGrid;
    op.evaluateOnGrid(testGeomDataOnGrid, trialGeomDataOnGrid, resultOnGrid);

    arma::Cube<ValueType> convertedResultOnGrid(1, 1, testPointCount * trialPointCount);
    for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint)
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
            convertedResultOnGrid(testPoint + trialPoint * testPointCount) =
                    resultOnGrid(0, testPoint, 0, trialPoint);

    // Collect data with evaluateAtPointPairs
    GeomData testGeomDataAtPointPairs, trialGeomDataAtPointPairs;
    testGeomDataAtPointPairs.globals.set_size(worldDim, testPointCount * trialPointCount);
    trialGeomDataAtPointPairs.globals.set_size(worldDim, testPointCount * trialPointCount);
    for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint)
        for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
            testGeomDataAtPointPairs.globals.col(testPoint + trialPoint * testPointCount) =
                    testGeomDataOnGrid.globals.col(testPoint);
            trialGeomDataAtPointPairs.globals.col(testPoint + trialPoint * testPointCount) =
                    trialGeomDataOnGrid.globals.col(trialPoint);
        }

    arma::Cube<ValueType> resultAtPointPairs;
    op.evaluateAtPointPairs(testGeomDataAtPointPairs, trialGeomDataAtPointPairs,
                            resultAtPointPairs);

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    resultAtPointPairs, convertedResultOnGrid, 1e-6));
}

BOOST_AUTO_TEST_SUITE_END()
