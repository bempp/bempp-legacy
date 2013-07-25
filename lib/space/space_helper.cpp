// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "space_helper.hpp"

#include "../common/acc.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::
getGlobalDofBoundingBoxes_defaultImplementation(
        const GridView& view,
        const std::vector<std::vector<LocalDof> >& global2localDofs,
        std::vector<BoundingBox<CoordinateType> >& bboxes)
{
    const IndexSet& indexSet = view.indexSet();
    const int elementCount = view.entityCount(0);

    std::vector<arma::Mat<CoordinateType> > elementCorners(elementCount);
    std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        const Geometry& geo = e.geometry();
        geo.getCorners(acc(elementCorners, index));
        it->next();
    }

    BoundingBox<CoordinateType> model;
    const CoordinateType maxCoord = std::numeric_limits<CoordinateType>::max();
    model.lbound.x = model.lbound.y = model.lbound.z = maxCoord;
    model.ubound.x = model.ubound.y = model.ubound.z = -maxCoord;

    const int globalDofCount_ = global2localDofs.size();
    bboxes.resize(globalDofCount_, model);
    for (int i = 0; i < globalDofCount_; ++i) {
        const std::vector<LocalDof>& localDofs = acc(global2localDofs, i);
        BoundingBox<CoordinateType>& bbox = acc(bboxes, i);
        for (int j = 0; j < localDofs.size(); ++j)
            extendBoundingBox(bbox, acc(elementCorners,
                                        acc(localDofs, j).entityIndex));
        assert(!localDofs.empty());
        setBoundingBoxReference<CoordinateType>(
                    bbox,
                    acc(elementCorners, localDofs[0].entityIndex).col(
                        localDofs[0].dofIndex));
    }

#ifndef NDEBUG
   for (size_t i = 0; i < globalDofCount_; ++i) {
       assert(bboxes[i].reference.x >= bboxes[i].lbound.x);
       assert(bboxes[i].reference.y >= bboxes[i].lbound.y);
       assert(bboxes[i].reference.z >= bboxes[i].lbound.z);
       assert(bboxes[i].reference.x <= bboxes[i].ubound.x);
       assert(bboxes[i].reference.y <= bboxes[i].ubound.y);
       assert(bboxes[i].reference.z <= bboxes[i].ubound.z);
   }
#endif // NDEBUG
}

template <typename BasisFunctionType>
void SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
        size_t flatLocalDofCount,
        const std::vector<std::vector<GlobalDofIndex> >& local2globalDofs,
        std::vector<LocalDof>& flatLocal2localDofs)
{
    flatLocal2localDofs.clear();
    flatLocal2localDofs.reserve(flatLocalDofCount);
    for (size_t e = 0; e < local2globalDofs.size(); ++e)
        for (size_t dof = 0; dof < acc(local2globalDofs, e).size(); ++dof)
            if (acc(acc(local2globalDofs, e), dof) >= 0)
                flatLocal2localDofs.push_back(LocalDof(e, dof));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(SpaceHelper);

} // namespace Bempp
