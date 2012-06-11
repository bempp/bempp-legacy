//// Copyright (C) 2011 by the BEM++ Authors
////
//// Permission is hereby granted, free of charge, to any person obtaining a copy
//// of this software and associated documentation files (the "Software"), to deal
//// in the Software without restriction, including without limitation the rights
//// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//// copies of the Software, and to permit persons to whom the Software is
//// furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included in
//// all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//// THE SOFTWARE.

//#include "fiber/laplace_3d_double_layer_potential_kernel.hpp"
//#include "fiber/laplace_3d_adjoint_double_layer_potential_kernel.hpp"
//#include "fiber/geometrical_data.hpp"
//#include "../type_template.hpp"
//#include "../check_arrays_are_close.hpp"

//#include <algorithm>
//#include <armadillo>
//#include <boost/test/unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/version.hpp>
//#include <complex>

//// Tests

//BOOST_AUTO_TEST_SUITE(Laplace3dAdjointDoubleLayerPotentialKernel)

//BOOST_AUTO_TEST_CASE_TEMPLATE(worldDimension_is_3, ValueType, kernel_types)
//{
//    Fiber::Laplace3dAdjointDoubleLayerPotentialKernel<ValueType> op;
//    BOOST_CHECK_EQUAL(op.worldDimension(), 3);
//}

//BOOST_AUTO_TEST_CASE_TEMPLATE(domainDimension_is_1, ValueType, kernel_types)
//{
//    Fiber::Laplace3dAdjointDoubleLayerPotentialKernel<ValueType> op;
//    BOOST_CHECK_EQUAL(op.domainDimension(), 1);
//}

//BOOST_AUTO_TEST_CASE_TEMPLATE(codomainDimension_is_1, ValueType, kernel_types)
//{
//    Fiber::Laplace3dAdjointDoubleLayerPotentialKernel<ValueType> op;
//    BOOST_CHECK_EQUAL(op.domainDimension(), 1);
//}

//BOOST_AUTO_TEST_CASE_TEMPLATE(addGeometricalDependencies_works,
//                              ValueType, kernel_types)
//{
//    Fiber::Laplace3dAdjointDoubleLayerPotentialKernel<ValueType> op;
//    int testGeomDeps = 1024, trialGeomDeps = 16; // random initial values
//    op.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

//    BOOST_CHECK(testGeomDeps & Fiber::GLOBALS);
//    BOOST_CHECK(testGeomDeps & Fiber::NORMALS);
//    BOOST_CHECK(trialGeomDeps & Fiber::GLOBALS);
//}

//BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_agrees_with_double_layer_potential,
//                              ValueType, kernel_types)
//{
//    // Check that DLP(x, y) = ADLP(y, x)

//    typedef Fiber::Laplace3dAdjointDoubleLayerPotentialKernel<ValueType> ADLPOperator;
//    ADLPOperator adlpOp;
//    Fiber::GeometricalData<typename ADLPOperator::CoordinateType> testGeomData, trialGeomData;
//    const int worldDim = 3;
//    const int testPointCount = 2, trialPointCount = 3;
//    testGeomData.globals.set_size(worldDim, testPointCount);
//    testGeomData.globals.fill(0.);
//    testGeomData.globals(0, 1) = 1.;

//    testGeomData.normals.set_size(worldDim, testPointCount);
//    testGeomData.normals.fill(0.);
//    testGeomData.normals(0, 0) = 1.;
//    testGeomData.normals(0, 1) = 1.;

//    trialGeomData.globals.set_size(worldDim, trialPointCount);
//    trialGeomData.globals.fill(1.);
//    trialGeomData.globals(0, 0) = 2.;
//    trialGeomData.globals(0, 1) = 3.;
//    trialGeomData.globals(0, 2) = 4.;

//    Fiber::_4dArray<ValueType> adlpResult;
//    adlpOp.evaluateOnGrid(testGeomData, trialGeomData, adlpResult);

//    typedef Fiber::Laplace3dDoubleLayerPotentialKernel<ValueType> DLPOperator;
//    DLPOperator dlpOp;
//    std::swap(testGeomData, trialGeomData);

//    Fiber::_4dArray<ValueType> dlpResult;
//    dlpOp.evaluateOnGrid(testGeomData, trialGeomData, dlpResult);

//    // swap test and trial point indexing
//    Fiber::_4dArray<ValueType> reorderedDlpResult(dlpResult.extent(2),
//                                                 dlpResult.extent(3),
//                                                 dlpResult.extent(0),
//                                                 dlpResult.extent(1));
//    for (int i0 = 0; i0 < dlpResult.extent(0); ++i0)
//        for (int i1 = 0; i1 < dlpResult.extent(1); ++i1)
//            for (int i2 = 0; i2 < dlpResult.extent(2); ++i2)
//                for (int i3 = 0; i3 < dlpResult.extent(3); ++i3)
//                    reorderedDlpResult(i2, i3, i0, i1) =
//                            dlpResult(i0, i1, i2, i3);

//    BOOST_CHECK(check_arrays_are_close<ValueType>(adlpResult, reorderedDlpResult, 1e-6));
//}

//BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_agrees_with_evaluateAtPointPairs,
//                              ValueType, kernel_types)
//{
//    typedef Fiber::Laplace3dAdjointDoubleLayerPotentialKernel<ValueType> Operator;
//    Operator op;

//    typedef Fiber::GeometricalData<typename Operator::CoordinateType> GeomData;

//    const int worldDim = 3;
//    const int testPointCount = 2, trialPointCount = 3;

//    // Collect data with evaluateOnGrid
//    GeomData testGeomDataOnGrid, trialGeomDataOnGrid;

//    testGeomDataOnGrid.globals.set_size(worldDim, testPointCount);
//    testGeomDataOnGrid.globals.fill(0.);
//    testGeomDataOnGrid.globals(0, 1) = 1.;

//    testGeomDataOnGrid.normals.set_size(worldDim, testPointCount);
//    testGeomDataOnGrid.normals.fill(0.);
//    testGeomDataOnGrid.normals(1, 0) = 1.;
//    testGeomDataOnGrid.normals(1, 1) = 1.;

//    trialGeomDataOnGrid.globals.set_size(worldDim, trialPointCount);
//    trialGeomDataOnGrid.globals.fill(1.);
//    trialGeomDataOnGrid.globals(0, 0) = 2.;
//    trialGeomDataOnGrid.globals(0, 1) = 3.;
//    trialGeomDataOnGrid.globals(0, 2) = 4.;

//    Fiber::_4dArray<ValueType> resultOnGrid;
//    op.evaluateOnGrid(testGeomDataOnGrid, trialGeomDataOnGrid, resultOnGrid);

//    arma::Cube<ValueType> convertedResultOnGrid(1, 1, testPointCount * trialPointCount);
//    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
//        for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
//            convertedResultOnGrid(testPoint + trialPoint * testPointCount) =
//                    resultOnGrid(0, testPoint, 0, trialPoint);

//    // Collect data with evaluateAtPointPairs
//    GeomData testGeomDataAtPointPairs, trialGeomDataAtPointPairs;
//    testGeomDataAtPointPairs.globals.set_size(worldDim, testPointCount * trialPointCount);
//    testGeomDataAtPointPairs.normals.set_size(worldDim, testPointCount * trialPointCount);
//    trialGeomDataAtPointPairs.globals.set_size(worldDim, testPointCount * trialPointCount);
//    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
//        for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
//            testGeomDataAtPointPairs.globals.col(testPoint + trialPoint * testPointCount) =
//                    testGeomDataOnGrid.globals.col(testPoint);
//            testGeomDataAtPointPairs.normals.col(testPoint + trialPoint * testPointCount) =
//                    testGeomDataOnGrid.normals.col(testPoint);
//            trialGeomDataAtPointPairs.globals.col(testPoint + trialPoint * testPointCount) =
//                    trialGeomDataOnGrid.globals.col(trialPoint);
//        }

//    arma::Cube<ValueType> resultAtPointPairs;
//    op.evaluateAtPointPairs(testGeomDataAtPointPairs, trialGeomDataAtPointPairs,
//                            resultAtPointPairs);

//    BOOST_CHECK(check_arrays_are_close<ValueType>(
//                    resultAtPointPairs, convertedResultOnGrid, 1e-6));
//}

//BOOST_AUTO_TEST_SUITE_END()
