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

#include "piecewise_constant_discontinuous_scalar_space_barycentric.hpp"

#include "space_helper.hpp"

#include "../common/acc.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_segment.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
PiecewiseConstantDiscontinuousScalarSpaceBarycentric(const shared_ptr<const Grid>& grid) :
    ScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
    m_segment(GridSegment::wholeGrid(*grid))
{
    assignDofsImpl(m_segment);
}

template <typename BasisFunctionType>
PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
PiecewiseConstantDiscontinuousScalarSpaceBarycentric(const shared_ptr<const Grid>& grid,
                             const GridSegment& segment) :
    ScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
    m_segment(segment)
{
    assignDofsImpl(m_segment);
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType> >
PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType> >& self) const
{
    if (self.get() != this)
        throw std::invalid_argument(
            "PiecewiseConstantDiscontinuousScalarSpaceBarycentric::discontinuousSpace(): "
            "argument should be a shared pointer to *this");
    return self;
}

template <typename BasisFunctionType>
bool
PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::isDiscontinuous() const
{
    return true;
}

template <typename BasisFunctionType>
int PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::domainDimension() const
{
    return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::codomainDimension() const
{
    return 1;
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType>&
PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::shapeset(
        const Entity<0>& element) const
{
    return m_shapeset;
}

template <typename BasisFunctionType>
ElementVariant PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.dim() == 1)
        return 2;
    if (type.isTriangle())
        return 3;
    else
        return 4;
}

template <typename BasisFunctionType>
bool PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
spaceIsCompatible(const Space<BasisFunctionType>& other) const
{

    if (other.grid().get()==this->grid().get()){
        return (other.spaceIdentifier()==this->spaceIdentifier());
    }
    else
        return false;

}


template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType> >
PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::barycentricSpace(
            const shared_ptr<const Space<BasisFunctionType> >& self) const {

    if (self.get()!=this)
        throw std::invalid_argument(
            "PiecewiseConstantDiscontinuousScalarSpaceBarycentric::barycentricSpace(): "
            "argument should be a shared pointer to *this");
     return self;
}


template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewiseConstantDiscontinuousScalarSpaceBarycentric::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::assignDofsImpl(
        const GridSegment& segment)
{

    const GridView& view = this->gridView();

    std::auto_ptr<GridView> viewCoarseGridPtr = this->grid()->levelView(0);
    const GridView& viewCoarseGrid = *viewCoarseGridPtr;

    const Mapper& elementMapper = view.elementMapper();
    const Mapper& elementMapperCoarseGrid = viewCoarseGrid.elementMapper();

    int elementCount = view.entityCount(0);
    int elementCountCoarseGrid = viewCoarseGrid.entityCount(0);

    // Assign gdofs to grid vertices (choosing only those that belong to
    // the selected grid segment)
    std::vector<int> continuousDofIndices(elementCountCoarseGrid, 0);
    segment.markExcludedEntities(0, continuousDofIndices);
    int continuousDofCount_ = 0;
    for (int elementIndex = 0; elementIndex < elementCountCoarseGrid; ++elementIndex)
        if (acc(continuousDofIndices, elementIndex) == 0) // not excluded
            acc(continuousDofIndices, elementIndex) = continuousDofCount_++;

    // (Re)initialise DOF maps
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_global2localDofs.clear();
    m_global2localDofs.reserve(elementCount);

    // Iterate over elements
    std::auto_ptr<EntityIterator<0> > itCoarseGrid = viewCoarseGrid.entityIterator<0>();
    int flatLocalDofCount_ = 0;
    while (!itCoarseGrid->finished()) {
        const Entity<0>& elementCoarseGrid = itCoarseGrid->entity();
        EntityIndex elementIndexCoarseGrid = elementMapperCoarseGrid.entityIndex(elementCoarseGrid);

        // Iterate through refined elements
        std::auto_ptr<EntityIterator<0> > sonIt = elementCoarseGrid.sonIterator(this->grid()->maxLevel());
        while (!sonIt->finished()){
            const Entity<0>& element = sonIt->entity();
            int elementIndex = elementMapper.entityIndex(element);
            std::vector<GlobalDofIndex>& globalDofs =
                    acc(m_local2globalDofs, elementIndex);
            int continuousDofIndex = acc(continuousDofIndices,elementIndexCoarseGrid);
            if (continuousDofIndex!= -1){
                globalDofs.push_back(flatLocalDofCount_);
                acc(m_global2localDofs,flatLocalDofCount_).push_back(
                            LocalDof(elementIndex, 0));
                ++flatLocalDofCount_;
            }
            else {
                globalDofs.push_back(-1);
            }
            sonIt->next();
        }
        itCoarseGrid->next();
    }

    // Initialize the container mapping the flat local dof indices to
    // local dof indices
    SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
                flatLocalDofCount_, m_local2globalDofs, m_flatLocal2localDofs);
}

template <typename BasisFunctionType>
size_t PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::flatLocalDofCount() const
{
    return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getGlobalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = this->gridView().elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
getGlobalDofInterpolationPoints(arma::Mat<CoordinateType>& points) const
{
    SpaceHelper<BasisFunctionType>::
            getGlobalDofInterpolationPoints_defaultImplementation(
                *this, points);
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::
getNormalsAtGlobalDofInterpolationPoints(arma::Mat<CoordinateType>& normals) const
{
    SpaceHelper<BasisFunctionType>::
            getNormalsAtGlobalDofInterpolationPoints_defaultImplementation(
                *this, normals);
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    std::vector<BoundingBox<CoordinateType> > bboxes;
    getGlobalDofBoundingBoxes(bboxes);

    positions.resize(bboxes.size());
    for (int i = 0; i < positions.size(); ++i)
        positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    std::vector<BoundingBox<CoordinateType> > bboxes;
    getFlatLocalDofBoundingBoxes(bboxes);

    positions.resize(bboxes.size());
    for (int i = 0; i < positions.size(); ++i)
        positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getGlobalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
    SpaceHelper<BasisFunctionType>::
            getGlobalDofBoundingBoxes_defaultImplementation(
                this->gridView(), m_global2localDofs, bboxes);
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getFlatLocalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
    getGlobalDofBoundingBoxes(bboxes);
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getGlobalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    SpaceHelper<BasisFunctionType>::
            getGlobalDofNormals_defaultImplementation(
                this->gridView(), m_global2localDofs, normals);
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    const int gridDim = this->domainDimension();
    const int worldDim = this->grid()->dimWorld();
    normals.resize(m_flatLocal2localDofs.size());

    const IndexSet& indexSet = this->gridView().indexSet();
    int elementCount = this->gridView().entityCount(0);

    arma::Mat<CoordinateType> elementNormals(worldDim, elementCount);
    std::auto_ptr<EntityIterator<0> > it = this->gridView().template entityIterator<0>();
    arma::Col<CoordinateType> center(gridDim);
    center.fill(0.5);
    arma::Col<CoordinateType> normal;
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        e.geometry().getNormals(center, normal);

        for (int dim = 0; dim < worldDim; ++dim)
            elementNormals(dim, index) = normal(dim);
        it->next();
    }

    if (gridDim == 1)
        for (size_t f = 0; f < m_flatLocal2localDofs.size(); ++f) {
            int elementIndex = m_flatLocal2localDofs[f].entityIndex;
            normals[f].x = elementNormals(0, elementIndex);
            normals[f].y = elementNormals(1, elementIndex);
            normals[f].z = 0.;
        }
    else // gridDim == 2
        for (size_t f = 0; f < m_flatLocal2localDofs.size(); ++f) {
            int elementIndex = m_flatLocal2localDofs[f].entityIndex;
            normals[f].x = elementNormals(0, elementIndex);
            normals[f].y = elementNormals(1, elementIndex);
            normals[f].z = elementNormals(2, elementIndex);
        }
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs) const
{
    dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseConstantDiscontinuousScalarSpaceBarycentric<BasisFunctionType>::dumpClusterIdsEx(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
        DofType dofType) const
{
    if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
        throw std::invalid_argument("PiecewiseConstantDiscontinuousScalarSpaceBarycentric::"
                                    "dumpClusterIds(): invalid DOF type");
    const size_t idCount = clusterIdsOfGlobalDofs.size();
    if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
            (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
        throw std::invalid_argument(
                "PiecewiseConstantDiscontinuousScalarSpaceBarycentric::dumpClusterIds(): "
                "clusterIds has incorrect length");

    std::auto_ptr<GridView> view = this->grid()->leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();
    arma::Row<double> data(idCount);
    for (size_t i = 0; i < idCount; ++i)
        data(i) = clusterIdsOfGlobalDofs[i];
    vtkWriter->addCellData(data, "ids");
    vtkWriter->write(fileName);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantDiscontinuousScalarSpaceBarycentric);

} // namespace Bempp


