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

#include "common_tests_for_spaces.hpp"
#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"

#include "assembly/assembly_options.hpp"
#include "assembly/context.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/l2_norm.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "common/scalar_traits.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "space/piecewise_polynomial_continuous_scalar_space.hpp"

#include <boost/test/floating_point_comparison.hpp>

using namespace Bempp;

template <typename ValueType_>
class LinearFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        result(0) = x + y + z;
    }
};

template <typename ValueType_>
class QuadraticFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        result(0) = x*x + y*y + z*z;
    }
};

template <typename ValueType_>
class CubicFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        result(0) = x*x*x + y*y*y + z*z*z;
    }
};

template <typename ValueType_>
class QuarticFunction
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }

    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        CoordinateType x = point(0), y = point(1), z = point(2);
        result(0) = x*x*x*x + y*y*y*y + z*z*z*z;
    }
};

// Tests

BOOST_AUTO_TEST_SUITE(PiecewisePolynomialContinuousScalarSpace_)

BOOST_AUTO_TEST_CASE_TEMPLATE(linear_function_can_be_expanded_in_linear_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 1));

    AccuracyOptions accuracyOptions;
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
                new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> function(
        context, space, space,
        surfaceNormalIndependentFunction(LinearFunction<RT>()));
    CT absoluteError, relativeError;
    estimateL2Error(
        function, surfaceNormalIndependentFunction(LinearFunction<RT>()),
        *quadStrategy, absoluteError, relativeError);
    BOOST_CHECK_SMALL(relativeError, 1000 * std::numeric_limits<CT>::epsilon() /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(quadratic_function_can_be_expanded_in_quadratic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2));

    AccuracyOptions accuracyOptions;
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> function(
        context, space, space,
        surfaceNormalIndependentFunction(QuadraticFunction<RT>()));
    CT absoluteError, relativeError;
    estimateL2Error(
        function, surfaceNormalIndependentFunction(QuadraticFunction<RT>()),
        *quadStrategy, absoluteError, relativeError);
    BOOST_CHECK_SMALL(relativeError, 1000 * std::numeric_limits<CT>::epsilon() /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(quadratic_function_can_be_expanded_in_cubic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3));

    AccuracyOptions accuracyOptions;
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> function(
        context, space, space,
        surfaceNormalIndependentFunction(QuadraticFunction<RT>()));
    CT absoluteError, relativeError;
    estimateL2Error(
        function, surfaceNormalIndependentFunction(QuadraticFunction<RT>()),
        *quadStrategy, absoluteError, relativeError);
    BOOST_CHECK_SMALL(relativeError, 1000 * std::numeric_limits<CT>::epsilon() /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(cubic_function_can_be_expanded_in_cubic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3));

    AccuracyOptions accuracyOptions;
    accuracyOptions.singleRegular.setRelativeQuadratureOrder(2);
    shared_ptr<NumericalQuadratureStrategy<BFT, RT> > quadStrategy(
        new NumericalQuadratureStrategy<BFT, RT>(accuracyOptions));
    AssemblyOptions assemblyOptions;
    assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
    shared_ptr<Context<BFT, RT> > context(
        new Context<BFT, RT>(quadStrategy, assemblyOptions));

    GridFunction<BFT, RT> function(
        context, space, space,
        surfaceNormalIndependentFunction(CubicFunction<RT>()));
    CT absoluteError, relativeError;
    estimateL2Error(
        function, surfaceNormalIndependentFunction(CubicFunction<RT>()),
        *quadStrategy, absoluteError, relativeError);
    BOOST_CHECK_SMALL(relativeError, 1000 * std::numeric_limits<CT>::epsilon() /* percent */);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_for_quadratic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2)));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_for_quadratic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2)));

    global2local_matches_local2global<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_for_cubic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3)));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_for_cubic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3)));

    global2local_matches_local2global<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_for_quadratic_space_and_segment, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2, segment)));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_for_cubic_space_and_segment, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3, segment)));

    global2local_matches_local2global<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_for_quadratic_space_and_segment_with_strictlyOnSegment_enabled, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(
             grid, 2, segment, true /* strictlyOnSegment */)));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_for_cubic_space_and_segment_with_strictlyOnSegment_enabled, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.2.msh", false /* verbose */);

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(
             grid, 3, segment, true /* strictlyOnSegment */)));

    global2local_matches_local2global<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(complement_is_really_a_complement_for_quadratic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2)));
    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space1(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2, segment)));
    GridSegment complement = segment.complement();
    shared_ptr<Space<BFT> > space2(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 2, complement)));

    complement_is_really_a_complement(space, space1, space2);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(complement_is_really_a_complement_for_cubic_space, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3)));
    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space1(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3, segment)));
    GridSegment complement = segment.complement();
    shared_ptr<Space<BFT> > space2(
        (new PiecewisePolynomialContinuousScalarSpace<BFT>(grid, 3, complement)));

    complement_is_really_a_complement(space, space1, space2);
}

BOOST_AUTO_TEST_SUITE_END()
