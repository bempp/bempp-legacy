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

#include "fiber/geometrical_data.hpp"
#include "fiber/laplace_3d_double_layer_potential_kernel_functor.hpp"
#include "fiber/standard_collection_of_kernels.hpp"

#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"

#include <algorithm>
#include "common/armadillo_fwd.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

BOOST_AUTO_TEST_SUITE(Laplace3dDoubleLayerPotentialKernelFunctor)

BOOST_AUTO_TEST_CASE_TEMPLATE(addGeometricalDependencies_works,
                              ValueType, kernel_types)
{
    Fiber::Laplace3dDoubleLayerPotentialKernelFunctor<ValueType> op;
    size_t testGeomDeps = 1024, trialGeomDeps = 16; // random initial values
    op.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

    BOOST_CHECK(testGeomDeps & Fiber::GLOBALS);
    BOOST_CHECK(trialGeomDeps & Fiber::GLOBALS);
    BOOST_CHECK(trialGeomDeps & Fiber::NORMALS);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_works_for_points_on_x_axis,
                              ValueType, kernel_types)
{
    typedef Fiber::Laplace3dDoubleLayerPotentialKernelFunctor<ValueType> Functor;
    typedef Fiber::StandardCollectionOfKernels<Functor> Kernels;
    Kernels kernels((Functor()));

    Fiber::GeometricalData<typename Functor::CoordinateType>
            testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 2, trialPointCount = 3;
    // Note: points lying on x axis only
    testGeomData.globals.set_size(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 1) = 1.;

    trialGeomData.globals.set_size(worldDim, trialPointCount);
    trialGeomData.globals.fill(0.);
    trialGeomData.globals(0, 0) = 2.;
    trialGeomData.globals(0, 1) = 3.;
    trialGeomData.globals(0, 2) = 4.;

    trialGeomData.normals.set_size(worldDim, trialPointCount);
    trialGeomData.normals.fill(0.);
    trialGeomData.normals(0, 0) = 1.;
    trialGeomData.normals(0, 1) = 1.;
    trialGeomData.normals(0, 2) = 1.;

    Fiber::CollectionOf4dArrays<ValueType> result;
    kernels.evaluateOnGrid(testGeomData, trialGeomData, result);

    Fiber::_4dArray<ValueType> expected(1, 1, testPointCount, trialPointCount);
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
        for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
            typename Functor::CoordinateType diff =
                    std::abs(testGeomData.globals(0, testPoint) -
                             trialGeomData.globals(0, trialPoint));
            expected(0, 0, testPoint, trialPoint) =
                    -(1. / (4. * M_PI)) / diff / diff;
        }

    BOOST_CHECK(check_arrays_are_close<ValueType>(result[0], expected, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_works_for_points_on_y_axis,
                              ValueType, kernel_types)
{
    typedef Fiber::Laplace3dDoubleLayerPotentialKernelFunctor<ValueType> Functor;
    typedef Fiber::StandardCollectionOfKernels<Functor> Kernels;
    Kernels kernels((Functor()));

    Fiber::GeometricalData<typename Functor::CoordinateType>
            testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 2, trialPointCount = 3;
    // Note: points lying on y axis only
    testGeomData.globals.set_size(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(1, 1) = 1.;

    trialGeomData.globals.set_size(worldDim, trialPointCount);
    trialGeomData.globals.fill(0.);
    trialGeomData.globals(1, 0) = 2.;
    trialGeomData.globals(1, 1) = 3.;
    trialGeomData.globals(1, 2) = 4.;

    trialGeomData.normals.set_size(worldDim, trialPointCount);
    trialGeomData.normals.fill(0.);
    trialGeomData.normals(1, 0) = 1.;
    trialGeomData.normals(1, 1) = 1.;
    trialGeomData.normals(1, 2) = 1.;

    Fiber::CollectionOf4dArrays<ValueType> result;
    kernels.evaluateOnGrid(testGeomData, trialGeomData, result);

    Fiber::_4dArray<ValueType> expected(1, 1, testPointCount, trialPointCount);
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
        for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
            typename Functor::CoordinateType diff =
                    std::abs(testGeomData.globals(1, testPoint) -
                             trialGeomData.globals(1, trialPoint));
            expected(0, 0, testPoint, trialPoint) =
                    -(1. / (4. * M_PI)) / diff / diff;
        }

    BOOST_CHECK(check_arrays_are_close<ValueType>(result[0], expected, 1e-6));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(evaluateOnGrid_agrees_with_evaluateAtPointPairs,
                              ValueType, kernel_types)
{
    typedef Fiber::Laplace3dDoubleLayerPotentialKernelFunctor<ValueType> Functor;
    typedef Fiber::StandardCollectionOfKernels<Functor> Kernels;
    Kernels kernels((Functor()));

    typedef Fiber::GeometricalData<typename Functor::CoordinateType> GeomData;

    const int worldDim = 3;
    const int testPointCount = 2, trialPointCount = 3;

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

    trialGeomDataOnGrid.normals.set_size(worldDim, trialPointCount);
    trialGeomDataOnGrid.normals.fill(0.);
    trialGeomDataOnGrid.normals(1, 0) = 1.;
    trialGeomDataOnGrid.normals(1, 1) = 1.;
    trialGeomDataOnGrid.normals(1, 2) = 1.;

    Fiber::CollectionOf4dArrays<ValueType> resultOnGrid;
    kernels.evaluateOnGrid(testGeomDataOnGrid, trialGeomDataOnGrid, resultOnGrid);

    arma::Col<ValueType> convertedResultOnGrid(testPointCount * trialPointCount);
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
        for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
            convertedResultOnGrid(testPoint + trialPoint * testPointCount) =
                    resultOnGrid[0](0, 0, testPoint, trialPoint);

    // Collect data with evaluateAtPointPairs
    GeomData testGeomDataAtPointPairs, trialGeomDataAtPointPairs;
    testGeomDataAtPointPairs.globals.set_size(worldDim, testPointCount * trialPointCount);
    trialGeomDataAtPointPairs.globals.set_size(worldDim, testPointCount * trialPointCount);
    trialGeomDataAtPointPairs.normals.set_size(worldDim, testPointCount * trialPointCount);
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
        for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
            testGeomDataAtPointPairs.globals.col(testPoint + trialPoint * testPointCount) =
                    testGeomDataOnGrid.globals.col(testPoint);
            trialGeomDataAtPointPairs.globals.col(testPoint + trialPoint * testPointCount) =
                    trialGeomDataOnGrid.globals.col(trialPoint);
            trialGeomDataAtPointPairs.normals.col(testPoint + trialPoint * testPointCount) =
                    trialGeomDataOnGrid.normals.col(trialPoint);
        }

    Fiber::CollectionOf3dArrays<ValueType> resultAtPointPairs;
    kernels.evaluateAtPointPairs(testGeomDataAtPointPairs, trialGeomDataAtPointPairs,
                                 resultAtPointPairs);
    arma::Col<ValueType> convertedResultAtPointPairs(testPointCount * trialPointCount);
    for (int point = 0; point < testPointCount * trialPointCount; ++point)
        convertedResultAtPointPairs(point) =
                resultAtPointPairs[0](0, 0, point);

    BOOST_CHECK(check_arrays_are_close<ValueType>(
                    convertedResultAtPointPairs, convertedResultOnGrid, 1e-6));
}

BOOST_AUTO_TEST_SUITE_END()
