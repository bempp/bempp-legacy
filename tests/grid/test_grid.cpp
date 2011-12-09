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

#include "test_grid.hpp"
#include "grid/grid_factory.hpp"
#include "grid/structured_grid_factory.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/version.hpp>

// Fixture member definitions

std::auto_ptr<Bempp::Grid> SimpleTriangularGridManager::createGrid()
{
    Bempp::GridParameters params;
    params.topology = Bempp::GridParameters::TRIANGULAR;

    const int dimGrid = 2;
    arma::Col<Bempp::ctype> lowerLeft(dimGrid);
    arma::Col<Bempp::ctype> upperRight(dimGrid);
    arma::Col<unsigned int> nElements(dimGrid);
    lowerLeft.fill(0);
    upperRight.fill(1);
    nElements(0) = N_ELEMENTS_X;
    nElements(1) = N_ELEMENTS_Y;

    return Bempp::GridFactory::createStructuredGrid(params, lowerLeft, upperRight, nElements);
}

std::auto_ptr<SimpleTriangularGridManager::DuneGrid> SimpleTriangularGridManager::createDuneGrid()
{
    const int dimGrid = 2;
    Dune::FieldVector<Bempp::ctype,dimGrid> duneLowerLeft;
    duneLowerLeft[0] = duneLowerLeft[1] = 0;
    Dune::FieldVector<Bempp::ctype,dimGrid> duneUpperRight;
    duneUpperRight[0] = duneUpperRight[1] = 1;
    Dune::array<unsigned int,dimGrid> duneNElements;
    duneNElements[0] = N_ELEMENTS_X;
    duneNElements[1] = N_ELEMENTS_Y;

    return Dune::BemppStructuredGridFactory<SimpleTriangularGridManager::DuneGrid>::
           createSimplexGrid(duneLowerLeft, duneUpperRight, duneNElements);
}

// Tests

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

// The test below fails because the numBoundarySegments method is not
// implemented in FoamGrid yet.

#if BOOST_VERSION >= 103500
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(boundarySegmentCount_agrees_with_Dune, 1)
#else
// See http://lists.boost.org/boost-users/2007/09/31144.php
#   ifdef _MSC_VER
#       pragma message("WARNING: Your version of Boost.Test does not register expected test failures correctly. Update to Boost 1.35 or newer.")
#   else // perhaps GCC
#       warning Your version of Boost.Test does not register expected test failures correctly. Update to Boost 1.35 or newer.
#   endif
#endif
BOOST_AUTO_TEST_CASE(boundarySegmentCount_agrees_with_Dune)
{
    BOOST_CHECK_EQUAL(bemppGrid->boundarySegmentCount(), duneGrid->numBoundarySegments());
}

BOOST_AUTO_TEST_SUITE_END()
