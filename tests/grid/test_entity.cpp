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
#include "grid/geometry.hpp"

#include <boost/test/unit_test.hpp>
#include "../num_template.hpp"

using namespace Bempp;

BOOST_FIXTURE_TEST_SUITE(Entity_Triangular, TriangularEntityManager)

// level()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(level_agrees_with_Dune_for_second_entity_on_level_0_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::Codim<codim>::EntityPointer duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    BOOST_CHECK_EQUAL((int)ep->entity().level(), (int)duneEp->level());
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(level_is_0_for_second_entity_on_level_0_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();

    BOOST_CHECK_EQUAL((int)ep->entity().level(), 0);
}

// type()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(type_agrees_with_Dune_for_second_entity_on_level_0_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;
    std::auto_ptr<EntityPointer<codim> > ep = getPointerToSecondEntityOnLevel0<codim>();
    typename DuneGrid::Codim<codim>::EntityPointer duneEp = getDunePointerToSecondEntityOnLevel0<codim>();

    BOOST_CHECK_EQUAL(ep->entity().type(), duneEp->type());
}

// subEntityCount()

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityCount_agrees_with_Dune_for_codimSub,
                                  T, list_1_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    int count = ep->entity().subEntityCount<codimSub>();

    typename DuneGrid::Codim<0>::EntityPointer duneEp = getDunePointerToSecondEntityOnLevel0<0>();
    int duneCount = duneEp->template count<codimSub>();

    BOOST_CHECK_EQUAL(count, duneCount);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityCount_is_zero_for_codimSub,
                                  T, list_0_and_3)
{
    const int codimSub = T::value;
    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    int count = ep->entity().subEntityCount<codimSub>();

    BOOST_CHECK_EQUAL(count, 0);
}

// subEntityIterator()

BOOST_AUTO_TEST_CASE(subEntityIterator_throws_for_codimSub_0)
{
    const int codimSub = 0;
    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    BOOST_CHECK_THROW(ep->entity().subEntityIterator<codimSub>(), std::logic_error);
}

BOOST_AUTO_TEST_CASE(subEntityIterator_does_not_throw_for_codimSub_1)
{
    const int codimSub = 1;
    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    BOOST_CHECK_NO_THROW(ep->entity().subEntityIterator<codimSub>());
}

BOOST_AUTO_TEST_CASE(subEntityIterator_does_not_throw_for_codimSub_2)
{
    const int codimSub = 2;
    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    BOOST_CHECK_NO_THROW(ep->entity().subEntityIterator<codimSub>());
}

BOOST_AUTO_TEST_CASE(subEntityIterator_throws_for_codimSub_3)
{
    const int codimSub = 3;
    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    BOOST_CHECK_THROW(ep->entity().subEntityIterator<codimSub>(), std::logic_error);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIterator_number_of_iterations_agrees_with_Dune_for_codimSub,
                                  T, list_1_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    std::auto_ptr<EntityIterator<codimSub> > it = ep->entity().subEntityIterator<codimSub>();
    int numIterations = 0;
    while (!it->finished()) {
        it->next();
        ++numIterations;
    }

    typename DuneGrid::Codim<0>::EntityPointer duneEp = getDunePointerToSecondEntityOnLevel0<0>();
    int duneNumIterations = duneEp->template count<codimSub>();

    BOOST_CHECK_EQUAL(numIterations, duneNumIterations);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIterator_number_of_iterations_is_3_for_codimSub,
                                  T, list_1_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
    std::auto_ptr<EntityIterator<codimSub> > it = ep->entity().subEntityIterator<codimSub>();
    int numIterations = 0;
    while (!it->finished()) {
        it->next();
        ++numIterations;
    }

    BOOST_CHECK_EQUAL(numIterations, 3);
}

/*
BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIterator_second_entity_agrees_with_Dune_for_codimSub,
                                  T, list_1_to_2)
{
    const size_t codimSub = T::value;

    arma::Col<double> elementCenter;
    {
        std::auto_ptr<EntityPointer<0> > ep = getPointerToSecondEntityOnLevel0<0>();
        std::auto_ptr<EntityIterator<codimSub> > it = ep->entity().subEntityIterator<codimSub>();
        it->next();
        const Entity<codimSub>& e = it->entity();
        const Geometry& geo = e.geometry();
        geo.getCenter(elementCenter);
    }

    Dune::FieldVector<DuneGrid::ctype, DuneGrid::dimensionworld> duneElementCenter;
    {
        typename DuneGrid::Codim<0>::EntityPointer duneEp = getDunePointerToSecondEntityOnLevel0<0>();
        typedef typename DuneGrid::Codim<codimSub> Codim;
        typename Codim::EntityPointer duneSubEp = duneEp->subEntity<codimSub>(1);
        const typename Codim::Entity& e = *duneSubEp;
        const typename Codim::Entity::Geometry& geo = e.geometry();
        duneElementCenter = geo.center();
    }

    BOOST_CHECK_EQUAL(elementCenter, duneElementCenter);
}
*/
BOOST_AUTO_TEST_SUITE_END()
