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

#include "../check_arrays_are_close.hpp"
#include "../type_template.hpp"
#include "../assembly/create_regular_grid.hpp"

#include "assembly/assembly_options.hpp"

#include "common/scalar_traits.hpp"

#include "grid/entity.hpp"
#include "grid/entity_iterator.hpp"
#include "grid/geometry.hpp"
#include "grid/grid.hpp"
#include "grid/grid_view.hpp"
#include "grid/grid_segment.hpp"
#include "grid/grid_factory.hpp"
#include "grid/index_set.hpp"
#include "grid/eigen_helpers.hpp"

#include <boost/test/floating_point_comparison.hpp>

using namespace Bempp;

namespace
{

template <int codim>
std::set<int> entitiesWithNonpositiveY(const GridView& view)
{
    std::set<int> result;
    std::unique_ptr<EntityIterator<codim> > it = view.entityIterator<codim>();
    const IndexSet& indexSet = view.indexSet();
    Vector<double> center;
    while (!it->finished()) {
        const Entity<codim>& entity = it->entity();
        entity.geometry().getCenter(Eigen::Ref<Vector<double>>(center));
        if (center(1) <= 0.)
            result.insert(indexSet.entityIndex(entity));
        it->next();
    }
    return result;
}

} // namespace

GridSegment gridSegmentWithPositiveY(const Grid& grid)
{
    boost::array<std::set<int>, 4> excludedEntities;
    std::unique_ptr<GridView> view = grid.leafView();
    excludedEntities[0] = entitiesWithNonpositiveY<0>(*view);
    excludedEntities[1] = entitiesWithNonpositiveY<1>(*view);
    excludedEntities[2] = entitiesWithNonpositiveY<2>(*view);

    return GridSegment(grid, excludedEntities[0], excludedEntities[1],
                       excludedEntities[2], excludedEntities[3]);
}

// Tests

BOOST_AUTO_TEST_SUITE(GridSegment_)

BOOST_AUTO_TEST_CASE(complement_works)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);

    GridSegment segment = gridSegmentWithPositiveX(*grid);
    GridSegment complement = segment.complement();

    std::unique_ptr<GridView> view = grid->leafView();
    for (int codim = 0; codim < 4; ++codim) {
        const int entityCount = view->entityCount(codim);
        for (int index = 0; index < entityCount; ++index) {
            BOOST_CHECK(segment.contains(codim, index) !=
                    complement.contains(codim, index));
        }
    }
}

BOOST_AUTO_TEST_CASE(union_works)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);

    GridSegment segmentA = gridSegmentWithPositiveX(*grid);
    GridSegment segmentB = gridSegmentWithPositiveY(*grid);
    GridSegment union_ = segmentA.union_(segmentB);

    std::unique_ptr<GridView> view = grid->leafView();
    for (int codim = 0; codim < 4; ++codim) {
        const int entityCount = view->entityCount(codim);
        for (int index = 0; index < entityCount; ++index) {
            bool oneOfSegmentsContainsIndex =
                    segmentA.contains(codim, index) ||
                    segmentB.contains(codim, index);
            bool unionContainsIndex = union_.contains(codim, index);
            BOOST_CHECK_EQUAL(oneOfSegmentsContainsIndex,
                              unionContainsIndex);
        }
    }
}

BOOST_AUTO_TEST_CASE(difference_works)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);

    GridSegment segmentA = gridSegmentWithPositiveX(*grid);
    GridSegment segmentB = gridSegmentWithPositiveY(*grid);
    GridSegment difference = segmentA.difference(segmentB);

    std::unique_ptr<GridView> view = grid->leafView();
    for (int codim = 0; codim < 4; ++codim) {
        const int entityCount = view->entityCount(codim);
        for (int index = 0; index < entityCount; ++index) {
            bool segmentAContainsIndexButSegmentBDoesNot =
                    segmentA.contains(codim, index) &&
                    !segmentB.contains(codim, index);
            bool differenceContainsIndex = difference.contains(codim, index);
            BOOST_CHECK_EQUAL(segmentAContainsIndexButSegmentBDoesNot,
                              differenceContainsIndex);
        }
    }
}

BOOST_AUTO_TEST_CASE(intersection_works)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "../../meshes/sphere-h-0.4.msh", false /* verbose */);

    GridSegment segmentA = gridSegmentWithPositiveX(*grid);
    GridSegment segmentB = gridSegmentWithPositiveY(*grid);
    GridSegment intersection = segmentA.intersection(segmentB);

    std::unique_ptr<GridView> view = grid->leafView();
    for (int codim = 0; codim < 4; ++codim) {
        const int entityCount = view->entityCount(codim);
        for (int index = 0; index < entityCount; ++index) {
            bool bothSegmentsContainIndex =
                    segmentA.contains(codim, index) &&
                    segmentB.contains(codim, index);
            bool intersectionContainsIndex = intersection.contains(codim, index);
            BOOST_CHECK_EQUAL(bothSegmentsContainIndex,
                              intersectionContainsIndex);
        }
    }
}

BOOST_AUTO_TEST_CASE(openDomain_works)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "meshes/cube-domains.msh", false /* verbose */);

    GridSegment segment = GridSegment::openDomain(*grid, 1);
    GridSegment complement = segment.complement();

    std::unique_ptr<GridView> view = grid->leafView();
    const IndexSet& indexSet = view->indexSet();
    std::unique_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        const int domain = e.domain();
        // Element
        const int index = indexSet.entityIndex(e);
        BOOST_CHECK_EQUAL(segment.contains(0, index), domain == 1);

        if (domain != 1) {
            // Check that none of the edges and vertices belongs to the segment
            const int edgeCount = e.subEntityCount<1>();
            for (int i = 0; i < edgeCount; ++i) {
                const int index = indexSet.subEntityIndex(e, i, 1);
                BOOST_CHECK(!segment.contains(1, index));
            }
            const int vertexCount = e.subEntityCount<2>();
            for (int i = 0; i < vertexCount; ++i) {
                const int index = indexSet.subEntityIndex(e, i, 2);
                BOOST_CHECK(!segment.contains(2, index));
            }
        }
        it->next();
    }
}

BOOST_AUTO_TEST_CASE(closedDomain_works)
{
    GridParameters params;
    params.topology = GridParameters::TRIANGULAR;
    shared_ptr<Grid> grid = GridFactory::importGmshGrid(
        params, "meshes/cube-domains.msh", false /* verbose */);

    GridSegment segment = GridSegment::closedDomain(*grid, 1);
    GridSegment complement = segment.complement();

    std::unique_ptr<GridView> view = grid->leafView();
    const IndexSet& indexSet = view->indexSet();
    std::unique_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        const int domain = e.domain();
        // Element
        const int index = indexSet.entityIndex(e);
        BOOST_CHECK_EQUAL(segment.contains(0, index), domain == 1);

        if (domain == 1) {
            // Check that all the edges and vertices belong to the segment
            const int edgeCount = e.subEntityCount<1>();
            for (int i = 0; i < edgeCount; ++i) {
                const int index = indexSet.subEntityIndex(e, i, 1);
                BOOST_CHECK(segment.contains(1, index));
            }
            const int vertexCount = e.subEntityCount<2>();
            for (int i = 0; i < vertexCount; ++i) {
                const int index = indexSet.subEntityIndex(e, i, 2);
                BOOST_CHECK(segment.contains(2, index));
            }
        }
        it->next();
    }
}

BOOST_AUTO_TEST_SUITE_END()
