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

#include "simple_triangular_grid_manager.hpp"
#include "grid/grid_factory.hpp"
#include "grid/structured_grid_factory.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>

BOOST_FIXTURE_TEST_SUITE(Grid, SimpleTriangularGridManager)

BOOST_AUTO_TEST_CASE(dimWorld_agrees_with_Dune)
{
    BOOST_CHECK_EQUAL(bemppGrid->dimWorld(), (int)duneGrid->dimensionworld);
}

BOOST_AUTO_TEST_CASE(dim_agrees_with_Dune)
{
    BOOST_CHECK_EQUAL(bemppGrid->dim(), (int)duneGrid->dimension);
}

BOOST_AUTO_TEST_CASE(maxLevel_agrees_with_Dune)
{
    BOOST_CHECK_EQUAL(bemppGrid->maxLevel(), duneGrid->maxLevel());
}

BOOST_AUTO_TEST_SUITE_END()
