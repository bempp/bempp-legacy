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
#include "grid/entity_iterator.hpp"
#include "grid/index_set.hpp"

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

BOOST_AUTO_TEST_SUITE_END()
