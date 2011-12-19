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

#include "test_index_set.hpp"
#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/index_set.hpp"

#include <stdexcept>
#include <boost/test/unit_test.hpp>
#include "../num_template.hpp"

using namespace Bempp;

BOOST_FIXTURE_TEST_SUITE(IndexSet_Triangular_Level0, TriangularLevel0IndexSetManager)

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIndex_agrees_with_Dune_for_second_entity_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    it->next();

    typename DuneGridView::Codim<codim>::Iterator duneIt = duneGridView.begin<codim>();
    ++duneIt;

    BOOST_CHECK_EQUAL(bemppIndexSet.entityIndex(it->entity()),
                      duneIndexSet.index(*duneIt));
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIndex_agrees_with_Dune_for_first_entity_of_codim_0_and_its_second_subentity_of_codim,
                                  T, list_0_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    typename DuneGridView::Codim<0>::Iterator duneIt = duneGridView.begin<0>();
    ++duneIt;

    bemppIndexSet.subEntityIndex(it->entity(), 0, codimSub);

    BOOST_CHECK_EQUAL(bemppIndexSet.subEntityIndex(it->entity(), 0, codimSub),
                      duneIndexSet.subIndex(*duneIt, 0, codimSub));
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIndex_and_entityIndex_agree_for_the_same_entity_of_codim,
                                  T, list_1_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    std::auto_ptr<EntityIterator<codimSub> > subIt = it->entity().subEntityIterator<codimSub>();
    subIt->next();

    int indexDirect = bemppIndexSet.entityIndex(subIt->entity());
    int indexIndirect = bemppIndexSet.subEntityIndex(it->entity(), 1, codimSub);

    BOOST_CHECK_EQUAL(indexDirect, indexIndirect);
}

BOOST_AUTO_TEST_CASE(subEntityIndex_returns_the_same_as_entityIndex_for_entity_of_codim_0)
{
    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    int indexDirect = bemppIndexSet.entityIndex(it->entity());
    int indexIndirect = bemppIndexSet.subEntityIndex(it->entity(), 0 /* i */, 0 /* codimSub */);

    BOOST_CHECK_EQUAL(indexDirect, indexIndirect);
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIndex_throws_for_invalid_codimension_and_codimSub,
                                  T, list_0_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    BOOST_CHECK_THROW(bemppIndexSet.subEntityIndex(it->entity(), 0, 3 /* codimSub */), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE(IndexSet_Triangular_Leaf, TriangularLeafIndexSetManager)

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(entityIndex_agrees_with_Dune_for_second_entity_of_codim,
                                  T, list_0_to_2)
{
    const int codim = T::value;

    std::auto_ptr<EntityIterator<codim> > it = bemppGridView->entityIterator<codim>();
    it->next();

    typename DuneGridView::Codim<codim>::Iterator duneIt = duneGridView.begin<codim>();
    ++duneIt;

    BOOST_CHECK_EQUAL(bemppIndexSet.entityIndex(it->entity()),
                      duneIndexSet.index(*duneIt));
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIndex_agrees_with_Dune_for_first_entity_of_codim_0_and_its_second_subentity_of_codim,
                                  T, list_0_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    typename DuneGridView::Codim<0>::Iterator duneIt = duneGridView.begin<0>();
    ++duneIt;

    bemppIndexSet.subEntityIndex(it->entity(), 0, codimSub);

    BOOST_CHECK_EQUAL(bemppIndexSet.subEntityIndex(it->entity(), 0, codimSub),
                      duneIndexSet.subIndex(*duneIt, 0, codimSub));
}

BOOST_AUTO_TEST_CASE_NUM_TEMPLATE(subEntityIndex_and_entityIndex_agree_for_the_same_entity_of_codim,
                                  T, list_1_to_2)
{
    const int codimSub = T::value;

    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    std::auto_ptr<EntityIterator<codimSub> > subIt = it->entity().subEntityIterator<codimSub>();
    subIt->next();

    int indexDirect = bemppIndexSet.entityIndex(subIt->entity());
    int indexIndirect = bemppIndexSet.subEntityIndex(it->entity(), 1, codimSub);

    BOOST_CHECK_EQUAL(indexDirect, indexIndirect);
}

BOOST_AUTO_TEST_CASE(subEntityIndex_returns_the_same_as_entityIndex_for_entity_of_codim_0)
{
    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    int indexDirect = bemppIndexSet.entityIndex(it->entity());
    int indexIndirect = bemppIndexSet.subEntityIndex(it->entity(), 0 /* i */, 0 /* codimSub */);

    BOOST_CHECK_EQUAL(indexDirect, indexIndirect);
}

#ifndef NDEBUG
BOOST_AUTO_TEST_CASE(subEntityIndex_throws_for_invalid_codimension)
{
    std::auto_ptr<EntityIterator<0> > it = bemppGridView->entityIterator<0>();
    it->next();

    BOOST_CHECK_THROW(bemppIndexSet.subEntityIndex(it->entity(), 0, 3 /* codimSub */), std::invalid_argument);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
