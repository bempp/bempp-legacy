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

#include "piecewise_linear_discontinuous_scalar_space.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_segment.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp
{

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
PiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid) :
    PiecewiseLinearScalarSpace<BasisFunctionType>(grid)
{
    GridSegment segment = GridSegment::wholeGrid(*grid);
    initialize(segment);
}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
PiecewiseLinearDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid,
                                        const GridSegment& segment) :
    PiecewiseLinearScalarSpace<BasisFunctionType>(grid)
{
    initialize(segment);
}

template <typename BasisFunctionType>
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
~PiecewiseLinearDiscontinuousScalarSpace()
{
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::initialize(
        const GridSegment& segment)
{
    const int gridDim = this->grid()->dim();
    if (gridDim != 1 && gridDim != 2)
        throw std::invalid_argument("PiecewiseLinearDiscontinuousScalarSpace::"
                                    "PiecewiseLinearDiscontinuousScalarSpace(): "
                                    "only 1- and 2-dimensional grids are supported");
    m_view = this->grid()->leafView();
    assignDofsImpl(segment);
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType> >
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType> >& self) const
{
    if (self.get() != this)
        throw std::invalid_argument(
            "PiecewiseLinearDiscontinuousScalarSpace::discontinuousSpace(): "
            "argument should be a shared pointer to *this");
    return self;
}

template <typename BasisFunctionType>
bool
PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::isDiscontinuous() const
{
    return true;
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::assignDofsImpl(
        const GridSegment& segment)
{
    // TODO: pay attention to the segment parameter
    const int gridDim = this->domainDimension();

    const Mapper& elementMapper = m_view->elementMapper();
    const IndexSet& indexSet = m_view->indexSet();

//    int globalDofCount_ = m_view->entityCount(this->grid()->dim());
    int elementCount = m_view->entityCount(0);

    // (Re)initialise DOF maps
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_global2localDofs.clear();
    // estimated number of global DOFs
    m_global2localDofs.reserve(4 * elementCount);
    // TODO: consider calling reserve(x) for each element of m_global2localDofs
    // with x being the typical number of elements adjacent to a vertex in a
    // grid of dimension gridDim

    // Iterate over elements
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    size_t globalDofCount = 0;
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        EntityIndex elementIndex = elementMapper.entityIndex(element);

        int vertexCount;
        if (gridDim == 1)
            vertexCount = element.template subEntityCount<1>();
        else // gridDim == 2
            vertexCount = element.template subEntityCount<2>();

        // List of global DOF indices corresponding to the local DOFs of the
        // current element
        std::vector<GlobalDofIndex>& globalDofs = m_local2globalDofs[elementIndex];
        globalDofs.reserve(vertexCount);
        for (int i = 0; i < vertexCount; ++i) {
            int vertexIndex = indexSet.subEntityIndex(element, i, gridDim);
            if (segment.contains(gridDim, vertexIndex)) {
                globalDofs.push_back(globalDofCount);
                std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, i));
                m_global2localDofs.push_back(localDofs);
                ++globalDofCount;
            } else
                globalDofs.push_back(-1);
        }
        it->next();
    }

    // Initialize the container mapping the flat local dof indices to
    // local dof indices
    m_flatLocal2localDofs.clear();
    m_flatLocal2localDofs.reserve(globalDofCount);
    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
        for (size_t dof = 0; dof < m_local2globalDofs[e].size(); ++dof)
            if (m_local2globalDofs[e][dof] >= 0)
                m_flatLocal2localDofs.push_back(LocalDof(e, dof));
}

template <typename BasisFunctionType>
size_t PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::flatLocalDofCount() const
{
    return globalDofCount();
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::getGlobalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i) {
        localDofs[i] = m_global2localDofs[globalDofs[i]];
        assert(localDofs[i].size() == 1);
    }
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    // This implementation assumes that the EntityIterator returns entities
    // ordered according to their indices
    const int worldDim = this->grid()->dimWorld();
    positions.resize(globalDofCount());

    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    arma::Mat<CoordinateType> corners;
    size_t globalDofIndex = 0;
    while (!it->finished())
    {
        const Entity<0>& e = it->entity();
        e.geometry().getCorners(corners);
        const size_t cornerCount = corners.n_cols;
        for (int corner = 0; corner < cornerCount; ++corner) {
            positions[globalDofIndex].x = corners(0, corner);
            positions[globalDofIndex].y = corners(1, corner);
            positions[globalDofIndex].z = (worldDim == 3) ? corners(2, corner) : 0.;
            ++globalDofIndex;
        }
        it->next();
    }
    assert(globalDofIndex == globalDofCount());
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    getGlobalDofPositions(positions);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
getGlobalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
   const int gridDim = this->domainDimension();
   const size_t globalDofCount_ = globalDofCount();
   bboxes.resize(globalDofCount_);

   BoundingBox<CoordinateType> model;
   model.lbound.x = std::numeric_limits<CoordinateType>::max();
   model.lbound.y = std::numeric_limits<CoordinateType>::max();
   model.lbound.z = std::numeric_limits<CoordinateType>::max();
   model.ubound.x = -std::numeric_limits<CoordinateType>::max();
   model.ubound.y = -std::numeric_limits<CoordinateType>::max();
   model.ubound.z = -std::numeric_limits<CoordinateType>::max();
   std::fill(bboxes.begin(), bboxes.end(), model);

   arma::Mat<CoordinateType> corners;

   if (gridDim != 2)
       throw std::runtime_error("PiecewiseLinearDiscontinuousScalarSpace::"
                                "getGlobalDofBoundingBoxes(): so far "
                                "implemented only for two-dimensional grids");

   std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
   size_t globalDofIndex = 0;
   while (!it->finished())
   {
       const Entity<0>& e = it->entity();
       const Geometry& geo = e.geometry();

       geo.getCorners(corners);
       const size_t cornerCount = corners.n_cols;
       for (size_t i = 0; i < cornerCount; ++i) {
           bboxes[globalDofIndex].reference.x = corners(0, i);
           bboxes[globalDofIndex].reference.y = corners(1, i);
           bboxes[globalDofIndex].reference.z = corners(2, i);
           for (size_t j = 0; j < cornerCount; ++j) {
               bboxes[globalDofIndex].lbound.x =
                   std::min(bboxes[globalDofIndex].lbound.x, corners(0, j));
               bboxes[globalDofIndex].lbound.y =
                   std::min(bboxes[globalDofIndex].lbound.y, corners(1, j));
               bboxes[globalDofIndex].lbound.z =
                   std::min(bboxes[globalDofIndex].lbound.z, corners(2, j));
               bboxes[globalDofIndex].ubound.x =
                   std::max(bboxes[globalDofIndex].ubound.x, corners(0, j));
               bboxes[globalDofIndex].ubound.y =
                   std::max(bboxes[globalDofIndex].ubound.y, corners(1, j));
               bboxes[globalDofIndex].ubound.z =
                   std::max(bboxes[globalDofIndex].ubound.z, corners(2, j));
           }
           ++globalDofIndex;
       }
       it->next();
   }
   assert(globalDofIndex == globalDofCount_);

#ifndef NDEBUG
   std::vector<Point3D<CoordinateType> > positions;
   getGlobalDofPositions(positions);
   for (size_t i = 0; i < globalDofCount_; ++i) {
       assert(bboxes[i].reference.x == positions[i].x);
       assert(bboxes[i].reference.y == positions[i].y);
       assert(bboxes[i].reference.z == positions[i].z);
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
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::
getFlatLocalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
    getGlobalDofBoundingBoxes(bboxes);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::getGlobalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    // This implementation assumes that the EntityIterator returns entities
    // ordered according to their indices
    const int gridDim = this->domainDimension();
    const int worldDim = this->grid()->dimWorld();
    normals.resize(globalDofCount());

    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    arma::Col<CoordinateType> center(gridDim);
    center.fill(0.5);
    arma::Col<CoordinateType> normal;

    size_t globalDofIndex = 0;
    while (!it->finished())
    {
        const Entity<0>& e = it->entity();
        e.geometry().getNormals(center, normal);
        int vertexCount;
        if (gridDim == 1)
            vertexCount = e.template subEntityCount<1>();
        else // gridDim == 2
            vertexCount = e.template subEntityCount<2>();
        for (int vertex = 0; vertex < vertexCount; ++vertex) {
            normals[globalDofIndex].x = normal(0);
            normals[globalDofIndex].y = normal(1);
            normals[globalDofIndex].z = (worldDim == 3) ? normal(2) : 0.;
            ++globalDofIndex;
        }
        it->next();
    }
    assert(globalDofIndex == globalDofCount());
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    getGlobalDofNormals(normals);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs) const
{
    dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>::dumpClusterIdsEx(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs,
        DofType dofType) const
{
    if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
        throw std::invalid_argument("PiecewiseLinearDiscontinuousScalarSpace::"
                                    "dumpClusterIds(): invalid DOF type");
    const size_t idCount = clusterIdsOfDofs.size();
    if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
            (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
        throw std::invalid_argument("PiecewiseLinearDiscontinuousScalarSpace::"
                                    "dumpClusterIds(): incorrect dimension");

    std::auto_ptr<GridView> view = this->grid()->leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();

    arma::Mat<double> data(idCount, globalDofCount());
    data.fill(0.);
    size_t row = 0;
    for (size_t id = 0; id < idCount; ++id) {
        bool exists = false;
        for (size_t fldof = 0; fldof < idCount; ++fldof) {
            if (clusterIdsOfDofs[fldof] == id) {
                LocalDof ldof = m_flatLocal2localDofs[fldof];
                GlobalDofIndex gdof = m_local2globalDofs[ldof.entityIndex][ldof.dofIndex];
                data(row, gdof) = 1;
                exists = true;
            }
        }
        if (!exists)
            data.shed_row(row); // very inefficient, of course
        else
            ++row;
    }
    std::cout << "about to write" <<std::endl;
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseLinearDiscontinuousScalarSpace);

} // namespace Bempp
