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

#include "test_entity.hpp"
#include "grid/armadillo_helpers.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "../num_template.hpp"

using namespace Bempp;

const double EPSILON = 1e-14;

BOOST_FIXTURE_TEST_SUITE(GridFactory_Triangular, TriangularEntityManager)

BOOST_AUTO_TEST_CASE(number_of_faces_is_correct)
{
    BOOST_CHECK_EQUAL(bemppGrid->levelView(0)->entityCount(0), 2 * N_ELEMENTS_X * N_ELEMENTS_Y);
}

BOOST_AUTO_TEST_CASE(number_of_vertices_is_correct)
{
    BOOST_CHECK_EQUAL(bemppGrid->levelView(0)->entityCount(2), (N_ELEMENTS_X + 1) * (N_ELEMENTS_Y + 1));
}

BOOST_AUTO_TEST_CASE(second_face_is_a_triangle)
{
    const int codim = 0;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    BOOST_CHECK(ep->entity().type().isTriangle());
}

BOOST_AUTO_TEST_CASE(elements_are_in_the_z_plane)
{
    std::auto_ptr<Bempp::GridView> bemppGridView = bemppGrid->levelView(0);
    std::auto_ptr<Bempp::EntityIterator<2> > it = bemppGridView->entityIterator<2>();

    ctype max_abs_z = 0.;
    arma::Col<ctype> center;
    while(!it->finished()) {
        it->entity().geometry().center(center);
        max_abs_z = std::max(max_abs_z, fabs(center(2)));
        it->next();
    }

    BOOST_CHECK_SMALL(max_abs_z, EPSILON);
}

BOOST_AUTO_TEST_CASE(elements_are_cover_the_unit_square)
{
    std::auto_ptr<Bempp::GridView> bemppGridView = bemppGrid->levelView(0);
    std::auto_ptr<Bempp::EntityIterator<2> > it = bemppGridView->entityIterator<2>();

    ctype min_x =  1e100;
    ctype max_x = -1e100;
    ctype min_y =  1e100;
    ctype max_y = -1e100;

    arma::Col<ctype> center;
    while(!it->finished()) {
        it->entity().geometry().center(center);
        max_x = std::max(max_x, center(0));
        max_y = std::max(max_y, center(1));
        min_x = std::min(min_x, center(0));
        min_y = std::min(min_y, center(1));
        it->next();
    }

    BOOST_CHECK_CLOSE(min_x, 0., EPSILON);
    BOOST_CHECK_CLOSE(max_x, 1., EPSILON);
    BOOST_CHECK_CLOSE(min_y, 0., EPSILON);
    BOOST_CHECK_CLOSE(max_y, 1., EPSILON);
}

BOOST_AUTO_TEST_CASE(jacobian_is_constant_everywhere_on_the_second_face)
{
    const int codim = 0;
    const int dimGlobal = DuneGrid::dimensionworld;
    const int dimLocal = DuneGrid::dimension - codim;
    const int nPoints = 5;

    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    const Geometry& geo = ep->entity().geometry();

    arma::Mat<ctype> local(dimLocal, nPoints);
    Dune::FieldVector<ctype, dimLocal> duneLocal;
    for (int j = 0; j < nPoints; ++j)
        for (int i = 0; i < dimLocal; ++i)
            local(i,j) = 0.1 * (i + 1) + 0.01 * (j + 1);

    arma::Row<ctype> intElement;
    geo.integrationElement(local, intElement);

    BOOST_CHECK_SMALL(stddev(intElement), EPSILON);
}

BOOST_AUTO_TEST_SUITE_END()
