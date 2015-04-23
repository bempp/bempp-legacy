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

#include "assembly/helmholtz_3d_operators_common.hpp"
#include "fiber/geometrical_data.hpp"
#include "fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp"
#include "fiber/modified_helmholtz_3d_single_layer_potential_kernel_interpolated_functor.hpp"
#include "fiber/default_collection_of_kernels.hpp"

#include "../type_template.hpp"
#include "../check_arrays_are_close.hpp"
#include "../random_arrays.hpp"

#include <algorithm>
#include "common/eigen_support.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>
#include <complex>

// Tests

BOOST_AUTO_TEST_SUITE(ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor)

BOOST_AUTO_TEST_CASE_TEMPLATE(agrees_with_noninterpolated_for_real_wave_number,
                              ValueType, kernel_types)
{
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>
            NoninterpolatedFunctor;
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<ValueType>
            InterpolatedFunctor;
    typedef Fiber::DefaultCollectionOfKernels<NoninterpolatedFunctor>
            NoninterpolatedKernels;
    typedef Fiber::DefaultCollectionOfKernels<InterpolatedFunctor>
            InterpolatedKernels;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    const ValueType waveNumber = 1;
    const CoordinateType wavelength = 2. * M_PI / std::abs(waveNumber);
    const double maxDist = 20.;
    const int interpPtsPerWavelength = Bempp::DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY;
    NoninterpolatedKernels noninterpKernels((NoninterpolatedFunctor(waveNumber)));
    InterpolatedKernels interpKernels((InterpolatedFunctor(waveNumber, maxDist,
                                                          interpPtsPerWavelength)));

    Fiber::GeometricalData<CoordinateType> testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 3, trialPointCount = 50;
    testGeomData.globals.resize(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 0) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(1, 1) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(2, 2) = wavelength / (2 * interpPtsPerWavelength);

    trialGeomData.globals = 0.5 * maxDist *
            generateRandomMatrix<CoordinateType>(worldDim, trialPointCount);
    trialGeomData.globals.col(0).fill(0.);
    trialGeomData.globals.block(0, 1, trialGeomData.globals.rows(), 10) *= 0.01; // to test well the area near origin
    trialGeomData.globals.block(0, 11, trialGeomData.globals.rows(), 10) *= 0.1;

    Fiber::CollectionOf4dArrays<ValueType> noninterpResult, interpResult;
    noninterpKernels.evaluateOnGrid(testGeomData, trialGeomData, noninterpResult);
    interpKernels.evaluateOnGrid(testGeomData, trialGeomData, interpResult);

    CoordinateType tol = 50 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(interpResult[0],
                                                  noninterpResult[0], tol));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(agrees_with_noninterpolated_for_negative_real_wave_number,
                              ValueType, kernel_types)
{
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>
            NoninterpolatedFunctor;
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<ValueType>
            InterpolatedFunctor;
    typedef Fiber::DefaultCollectionOfKernels<NoninterpolatedFunctor>
            NoninterpolatedKernels;
    typedef Fiber::DefaultCollectionOfKernels<InterpolatedFunctor>
            InterpolatedKernels;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    const ValueType waveNumber = -1;
    const CoordinateType wavelength = 2. * M_PI / std::abs(waveNumber);
    const double maxDist = 20.;
    const int interpPtsPerWavelength = Bempp::DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY;
    NoninterpolatedKernels noninterpKernels((NoninterpolatedFunctor(waveNumber)));
    InterpolatedKernels interpKernels((InterpolatedFunctor(waveNumber, maxDist,
                                                          interpPtsPerWavelength)));

    Fiber::GeometricalData<CoordinateType> testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 3, trialPointCount = 50;
    testGeomData.globals.resize(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 0) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(1, 1) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(2, 2) = wavelength / (2 * interpPtsPerWavelength);

    trialGeomData.globals = 0.5 * maxDist *
            generateRandomMatrix<CoordinateType>(worldDim, trialPointCount);
    trialGeomData.globals.col(0).fill(0.);
    trialGeomData.globals.block(0, 1, trialGeomData.globals.rows(), 10) *= 0.01; // to test well the area near origin
    trialGeomData.globals.block(0, 11, trialGeomData.globals.rows(), 10) *= 0.1;

    Fiber::CollectionOf4dArrays<ValueType> noninterpResult, interpResult;
    noninterpKernels.evaluateOnGrid(testGeomData, trialGeomData, noninterpResult);
    interpKernels.evaluateOnGrid(testGeomData, trialGeomData, interpResult);

    CoordinateType tol = 50 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(interpResult[0],
                                                  noninterpResult[0], tol));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(agrees_with_noninterpolated_for_negative_imag_wave_number,
                              ValueType, complex_kernel_types)
{
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>
            NoninterpolatedFunctor;
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<ValueType>
            InterpolatedFunctor;
    typedef Fiber::DefaultCollectionOfKernels<NoninterpolatedFunctor>
            NoninterpolatedKernels;
    typedef Fiber::DefaultCollectionOfKernels<InterpolatedFunctor>
            InterpolatedKernels;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    const ValueType waveNumber(0., -1.);
    const CoordinateType wavelength = 2. * M_PI / std::abs(waveNumber);
    const double maxDist = 20.;
    const int interpPtsPerWavelength = Bempp::DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY;
    NoninterpolatedKernels noninterpKernels((NoninterpolatedFunctor(waveNumber)));
    InterpolatedKernels interpKernels((InterpolatedFunctor(waveNumber, maxDist,
                                                          interpPtsPerWavelength)));

    Fiber::GeometricalData<CoordinateType> testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 3, trialPointCount = 30;
    testGeomData.globals.resize(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 0) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(1, 1) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(2, 2) = wavelength / (2 * interpPtsPerWavelength);

    trialGeomData.globals = 0.5 * maxDist *
            generateRandomMatrix<CoordinateType>(worldDim, trialPointCount);
    trialGeomData.globals.col(0).fill(0.);
    trialGeomData.globals.block(0, 1, trialGeomData.globals.rows(), 10) *= 0.01; // to test well the area near origin
    trialGeomData.globals.block(0, 11, trialGeomData.globals.rows(), 10) *= 0.1;

    Fiber::CollectionOf4dArrays<ValueType> noninterpResult, interpResult;
    noninterpKernels.evaluateOnGrid(testGeomData, trialGeomData, noninterpResult);
    interpKernels.evaluateOnGrid(testGeomData, trialGeomData, interpResult);

    CoordinateType tol = 50 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(interpResult[0],
                                                  noninterpResult[0], tol));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(agrees_with_noninterpolated_for_imag_wave_number,
                              ValueType, complex_kernel_types)
{
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>
            NoninterpolatedFunctor;
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<ValueType>
            InterpolatedFunctor;
    typedef Fiber::DefaultCollectionOfKernels<NoninterpolatedFunctor>
            NoninterpolatedKernels;
    typedef Fiber::DefaultCollectionOfKernels<InterpolatedFunctor>
            InterpolatedKernels;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    const ValueType waveNumber(0., 1.);
    const CoordinateType wavelength = 2. * M_PI / std::abs(waveNumber);
    const double maxDist = 20.;
    const int interpPtsPerWavelength = Bempp::DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY;
    NoninterpolatedKernels noninterpKernels((NoninterpolatedFunctor(waveNumber)));
    InterpolatedKernels interpKernels((InterpolatedFunctor(waveNumber, maxDist,
                                                          interpPtsPerWavelength)));

    Fiber::GeometricalData<CoordinateType> testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 3, trialPointCount = 30;
    testGeomData.globals.resize(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 0) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(1, 1) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(2, 2) = wavelength / (2 * interpPtsPerWavelength);

    trialGeomData.globals = 0.5 * maxDist *
            generateRandomMatrix<CoordinateType>(worldDim, trialPointCount);
    trialGeomData.globals.col(0).fill(0.);
    trialGeomData.globals.block(0, 1, trialGeomData.globals.rows(), 10) *= 0.01; // to test well the area near origin
    trialGeomData.globals.block(0, 11, trialGeomData.globals.rows(), 10) *= 0.1;

    Fiber::CollectionOf4dArrays<ValueType> noninterpResult, interpResult;
    noninterpKernels.evaluateOnGrid(testGeomData, trialGeomData, noninterpResult);
    interpKernels.evaluateOnGrid(testGeomData, trialGeomData, interpResult);

    CoordinateType tol = 50 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(interpResult[0],
                                                  noninterpResult[0], tol));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(agrees_with_noninterpolated_for_complex_wave_number,
                              ValueType, complex_kernel_types)
{
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>
            NoninterpolatedFunctor;
    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<ValueType>
            InterpolatedFunctor;
    typedef Fiber::DefaultCollectionOfKernels<NoninterpolatedFunctor>
            NoninterpolatedKernels;
    typedef Fiber::DefaultCollectionOfKernels<InterpolatedFunctor>
            InterpolatedKernels;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    const ValueType waveNumber(1., 1.);
    const CoordinateType wavelength = 2. * M_PI / std::abs(waveNumber);
    const double maxDist = 20.;
    const int interpPtsPerWavelength = Bempp::DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY;
    NoninterpolatedKernels noninterpKernels((NoninterpolatedFunctor(waveNumber)));
    InterpolatedKernels interpKernels((InterpolatedFunctor(waveNumber, maxDist,
                                                          interpPtsPerWavelength)));

    Fiber::GeometricalData<CoordinateType> testGeomData, trialGeomData;
    const int worldDim = 3;
    const int testPointCount = 3, trialPointCount = 30;
    testGeomData.globals.resize(worldDim, testPointCount);
    testGeomData.globals.fill(0.);
    testGeomData.globals(0, 0) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(1, 1) = wavelength / (2 * interpPtsPerWavelength);
    testGeomData.globals(2, 2) = wavelength / (2 * interpPtsPerWavelength);

    trialGeomData.globals = 0.5 * maxDist *
            generateRandomMatrix<CoordinateType>(worldDim, trialPointCount);
    trialGeomData.globals.col(0).fill(0.);
    trialGeomData.globals.block(0, 1, trialGeomData.globals.rows(), 10) *= 0.01; // to test well the area near origin
    trialGeomData.globals.block(0, 11, trialGeomData.globals.rows(), 10) *= 0.1;

    Fiber::CollectionOf4dArrays<ValueType> noninterpResult, interpResult;
    noninterpKernels.evaluateOnGrid(testGeomData, trialGeomData, noninterpResult);
    interpKernels.evaluateOnGrid(testGeomData, trialGeomData, interpResult);

    CoordinateType tol = 50 * std::numeric_limits<CoordinateType>::epsilon();
    BOOST_CHECK(check_arrays_are_close<ValueType>(interpResult[0],
                                                  noninterpResult[0], tol));
}

BOOST_AUTO_TEST_SUITE_END()
