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

#include "test_grid_view.hpp"
#include "grid/armadillo_helpers.hpp"
#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#include "../num_template.hpp"

using namespace Bempp;

BOOST_FIXTURE_TEST_SUITE(GridView_Triangular_Level0, TriangularLevel0GridViewManager)

// entityCount()

BOOST_AUTO_TEST_CASE(entityCount_is_zero_for_codim_3)
{
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(3), 0u);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityCount_agrees_with_Dune_for_codim,
                                  T, list_0_to_3)
{
    const int codim = T::value;
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(codim), (size_t) duneGridView.size(codim));
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityCount_agrees_with_Dune_for_simplex_of_dim,
                                  T, list_0_to_3)
{
    const int dim = T::value;
    const GeometryType type(GeometryType::simplex, dim);
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(type), (size_t) duneGridView.size(type));
}

BOOST_AUTO_TEST_CASE(entityCount_agrees_with_Dune_for_cube_of_dim_2)
{
    const GeometryType type(GeometryType::cube, 2);
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(type), (size_t) duneGridView.size(type));
}

BOOST_AUTO_TEST_CASE(entityCount_is_zero_for_cube_of_dim_2)
{
    const GeometryType type(GeometryType::cube, 2);
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(type), 0u);
}

// entityIterator()

BOOST_AUTO_TEST_CASE(entityIterator_throws_for_codim_3)
{
    BOOST_CHECK_THROW(bemppGridView->entityIterator<3>(), std::logic_error);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_does_not_throw_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    BOOST_CHECK_NO_THROW(bemppGridView->entityIterator<codim>());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_number_of_iterations_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    int numIterations = 0;
    while (!it->finished()) {
        it->next();
        ++numIterations;
    }

    typedef typename DuneGridView::Codim<codim>::Iterator DuneIterator;
    DuneIterator duneIt = duneGridView.begin<codim>();
    const DuneIterator end = duneGridView.end<codim>();
    int duneNumIterations = 0;
    while (duneIt != end) {
        ++duneIt;
        ++duneNumIterations;
    }

    BOOST_CHECK_EQUAL(numIterations, duneNumIterations);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_finished_does_not_return_true_immediately_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    BOOST_CHECK(!it->finished());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_second_entity_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    arma::Col<double> elementCenter;
    {
        std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
        it->next();
        const Entity<codim>& e = it->entity();
        const Geometry& geo = e.geometry();
        geo.getCenter(elementCenter);
    }

    Dune::FieldVector<DuneGrid::ctype, DuneGrid::dimensionworld> duneElementCenter;
    {
        typedef typename DuneGridView::Codim<codim> Codim;
        typename Codim::Iterator duneIt = duneGridView.begin<codim>();
        ++duneIt;
        const typename Codim::Entity& e = *duneIt;
        const typename Codim::Entity::Geometry& geo = e.geometry();
        duneElementCenter = geo.center();
    }

    BOOST_CHECK_EQUAL(elementCenter, duneElementCenter);
}

// containsEntity()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(containsEntity_returns_true_for_second_entity_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    it->next();
    BOOST_CHECK(bemppGridView->containsEntity(it->entity()));
}

BOOST_AUTO_TEST_SUITE_END()





BOOST_FIXTURE_TEST_SUITE(GridView_Triangular_Leaf, TriangularLeafGridViewManager)

// Contents of this test suite are identical with those of the previous one -- the only difference is the fixture.
// I can't see any elegant way of avoiding this duplication, though.

// entityCount()

BOOST_AUTO_TEST_CASE(entityCount_is_zero_for_codim_3)
{
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(3), 0u);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityCount_agrees_with_Dune_for_codim,
                                  T, list_0_to_3)
{
    const int codim = T::value;
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(codim), (size_t) duneGridView.size(codim));
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityCount_agrees_with_Dune_for_simplex_of_dim,
                                  T, list_0_to_3)
{
    const int dim = T::value;
    const GeometryType type(GeometryType::simplex, dim);
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(type), (size_t) duneGridView.size(type));
}

BOOST_AUTO_TEST_CASE(entityCount_agrees_with_Dune_for_cube_of_dim_2)
{
    const GeometryType type(GeometryType::cube, 2);
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(type), (size_t) duneGridView.size(type));
}

// This test fails -- there is a bug in FoamGridLeafIndexSet::size(GeometryType)
// (it just checks dimensions, the geometry type (simplex/cube/...) is not checked)
#if BOOST_VERSION >= 103500
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(entityCount_is_zero_for_cube_of_dim_2, 1)
#else
// See http://lists.boost.org/boost-users/2007/09/31144.php
#   ifdef _MSC_VER
#       pragma message("WARNING: Your version of Boost.Test does not register expected test failures correctly. Update to Boost 1.35 or newer.")
#   else // perhaps GCC
#       warning Your version of Boost.Test does not register expected test failures correctly. Update to Boost 1.35 or newer.
#   endif
#endif
BOOST_AUTO_TEST_CASE(entityCount_is_zero_for_cube_of_dim_2)
{
    const GeometryType type(GeometryType::cube, 2);
    BOOST_CHECK_EQUAL(bemppGridView->entityCount(type), 0u);
}

// entityIterator()

BOOST_AUTO_TEST_CASE(entityIterator_throws_for_codim_3)
{
    BOOST_CHECK_THROW(bemppGridView->entityIterator<3>(), std::logic_error);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_does_not_throw_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    BOOST_CHECK_NO_THROW(bemppGridView->entityIterator<codim>());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_number_of_iterations_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    int numIterations = 0;
    while (!it->finished()) {
        it->next();
        ++numIterations;
    }

    typedef typename DuneGridView::Codim<codim>::Iterator DuneIterator;
    DuneIterator duneIt = duneGridView.begin<codim>();
    const DuneIterator end = duneGridView.end<codim>();
    int duneNumIterations = 0;
    while (duneIt != end) {
        ++duneIt;
        ++duneNumIterations;
    }

    BOOST_CHECK_EQUAL(numIterations, duneNumIterations);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_finished_does_not_return_true_immediately_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    // Testing behaviour at failure
    //    if (codim == 1)
    //        BOOST_CHECK(0);
    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    BOOST_CHECK(!it->finished());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIterator_second_entity_agrees_with_Dune_for_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    arma::Col<double> elementCenter;
    {
        std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
        it->next();
        const Entity<codim>& e = it->entity();
        const Geometry& geo = e.geometry();
        geo.getCenter(elementCenter);
    }

    Dune::FieldVector<DuneGrid::ctype, DuneGrid::dimensionworld> duneElementCenter;
    {
        typedef typename DuneGridView::Codim<codim> Codim;
        typename Codim::Iterator duneIt = duneGridView.begin<codim>();
        ++duneIt;
        const typename Codim::Entity& e = *duneIt;
        const typename Codim::Entity::Geometry& geo = e.geometry();
        duneElementCenter = geo.center();
    }

    BOOST_CHECK_EQUAL(elementCenter, duneElementCenter);
}

// containsEntity()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(containsEntity_returns_true_for_second_entity_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    it->next();
    BOOST_CHECK(bemppGridView->containsEntity(it->entity()));
}

BOOST_AUTO_TEST_SUITE_END()


// IdSet: entityId should return the same value for dune and bempp for e.g. third entity with various codims. It should throw for codim 3.
// IndexSet: the same thing.

// Iterator: should iterate the same number of times.
// On start should return the same entity as Dune's (e.g. test by comparing centres)
// After a fixed number of iterations should also return the same entity.
// What should it do for an empty grid? (throw some exception probably...)
// StructuredGridFactory should return a grid with the correct number of elements.
// Geometry should return correct values of jacobians etc. for both one and multiple requested points.
