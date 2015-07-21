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
#include "../assembly/create_regular_grid.hpp"

#include "common/scalar_traits.hpp"

#include "grid/grid.hpp"
#include "grid/grid_factory.hpp"

#include "space/piecewise_linear_continuous_scalar_space.hpp"

#include <boost/type_traits/is_complex.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace Bempp;

// Tests

BOOST_AUTO_TEST_SUITE(PiecewiseLinearContinuousScalarSpace_)

BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.1.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid)));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.1.msh", false /* verbose */);

    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid)));

    global2local_matches_local2global<BFT>(*space);
}


/*
BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_for_segment, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.1.msh", false );

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid, segment)));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_for_segment, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.1.msh", false );

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid, segment)));

    global2local_matches_local2global<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(local2global_matches_global2local_for_segment_with_strictlyOnSegment_enabled, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.1.msh", false );

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(
             grid, segment, true )));

    local2global_matches_global2local<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(global2local_matches_local2global_for_segment_with_strictlyOnSegment_enabled, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.1.msh", false );

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(
             grid, segment, true )));

    global2local_matches_local2global<BFT>(*space);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(complement_is_really_a_complement_, ResultType, result_types)
{
    typedef ResultType RT;
    typedef typename ScalarTraits<RT>::RealType BFT;
    typedef typename ScalarTraits<RT>::RealType CT;

    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false );

    shared_ptr<Space<BFT> > space(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid)));
    GridSegment segment = gridSegmentWithPositiveX(*grid);
    shared_ptr<Space<BFT> > space1(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid, segment)));
    GridSegment complement = segment.complement();
    shared_ptr<Space<BFT> > space2(
        (new PiecewiseLinearContinuousScalarSpace<BFT>(grid, complement)));

    complement_is_really_a_complement(space, space1, space2);
}
*/
BOOST_AUTO_TEST_SUITE_END()
