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

#include "piecewise_linear_continuous_scalar_space.hpp"

#include "piecewise_linear_discontinuous_scalar_space.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp
{

template <typename BasisFunctionType>
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
PiecewiseLinearContinuousScalarSpace(const shared_ptr<const Grid>& grid) :
    PiecewiseLinearScalarSpace<BasisFunctionType>(grid), m_flatLocalDofCount(0)
{
    const int gridDim = grid->dim();
    if (gridDim != 1 && gridDim != 2)
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "PiecewiseLinearContinuousScalarSpace(): "
                                    "only 1- and 2-dimensional grids are supported");
    m_view = grid->leafView();
    assignDofsImpl();
}

template <typename BasisFunctionType>
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
~PiecewiseLinearContinuousScalarSpace()
{
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType> >
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType> >& self) const
{
    if (!m_discontinuousSpace) {
        tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
        typedef PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType>
                DiscontinuousSpace;
        if (!m_discontinuousSpace)
            m_discontinuousSpace.reset(new DiscontinuousSpace(this->grid()));
    }
    return m_discontinuousSpace;
}

template <typename BasisFunctionType>
bool
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::isDiscontinuous() const
{
    return false;
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::assignDofsImpl()
{
    const int gridDim = this->domainDimension();

    const Mapper& elementMapper = m_view->elementMapper();
    const IndexSet& indexSet = m_view->indexSet();

    // Global DOF numbers will be identical with vertex indices.
    // Thus, the will be as many global DOFs as there are vertices.
    int globalDofCount_ = m_view->entityCount(this->grid()->dim());
    int elementCount = m_view->entityCount(0);

//    // DEBUG
//    {
//        std::cout << "Vertices:\n" << std::endl;
//        std::auto_ptr<EntityIterator<2> > vit = m_view->entityIterator<2>();
//        const IndexSet& indexSet = m_view->indexSet();
//        while (!vit->finished())
//        {
//            arma::Col<BasisFunctionType> vertex;
//            vit->entity().geometry().center(vertex);
//            std::cout << indexSet.entityIndex(vit->entity()) << ": "
//                      << vertex(0) << " "
//                      << vertex(1) << " "
//                      << vertex(2) << std::endl;
//            vit->next();
//        }
//    }

    // (Re)initialise DOF maps
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_global2localDofs.clear();
    m_global2localDofs.resize(globalDofCount_);
    // TODO: consider calling reserve(x) for each element of m_global2localDofs
    // with x being the typical number of elements adjacent to a vertex in a
    // grid of dimension gridDim

    // Iterate over elements
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    m_flatLocalDofCount = 0;
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        EntityIndex elementIndex = elementMapper.entityIndex(element);

        int vertexCount;
        if (gridDim == 1)
            vertexCount = element.template subEntityCount<1>();
        else // gridDim == 2
            vertexCount = element.template subEntityCount<2>();
        m_flatLocalDofCount += vertexCount;

        // List of global DOF indices corresponding to the local DOFs of the
        // current element
        std::vector<GlobalDofIndex>& globalDofs = m_local2globalDofs[elementIndex];
        globalDofs.resize(vertexCount);
        for (int i = 0; i < vertexCount; ++i)
        {
            GlobalDofIndex globalDofIndex = indexSet.subEntityIndex(element, i, gridDim);
            globalDofs[i] = globalDofIndex;
            m_global2localDofs[globalDofIndex].push_back(LocalDof(elementIndex, i));
        }
        it->next();
    }

    // Initialize the container mapping the flat local dof indices to
    // local dof indices
    m_flatLocal2localDofs.clear();
    m_flatLocal2localDofs.reserve(m_flatLocalDofCount);
    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
        for (size_t dof = 0; dof < m_local2globalDofs[e].size(); ++dof)
            m_flatLocal2localDofs.push_back(LocalDof(e, dof));
}

template <typename BasisFunctionType>
size_t PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::flatLocalDofCount() const
{
// This is the correct implementation. Include it once the
// bug in FoamGrid is fixes.
//    if (gridDim == 1)
//        return m_view->entityCount(0) * 2;
//    else // gridDim == 2
//        return m_view->entityCount(GeometryType(GeometryType::cube, 2)) * 4 +
//                m_view->entityCount(GeometryType(GeometryType::simplex, 2)) * 3;

    return m_flatLocalDofCount;
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getGlobalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = this->domainDimension();
    const int globalDofCount_ = globalDofCount();
    positions.resize(globalDofCount_);

    const IndexSet& indexSet = m_view->indexSet();

    if (gridDim == 1)
    {
        std::auto_ptr<EntityIterator<1> > it = m_view->entityIterator<1>();
        while (!it->finished())
        {
            const Entity<1>& e = it->entity();
            int index = indexSet.entityIndex(e);
            arma::Col<CoordinateType> vertex;
            e.geometry().getCenter(vertex);

            positions[index].x = vertex(0);
            positions[index].y = vertex(1);
            positions[index].z = 0.;
            it->next();
        }
    }
    else // gridDim == 2
    {
        std::auto_ptr<EntityIterator<2> > it = m_view->entityIterator<2>();
        while (!it->finished())
        {
            const Entity<2>& e = it->entity();
            int index = indexSet.entityIndex(e);
            arma::Col<CoordinateType> vertex;
            e.geometry().getCenter(vertex);

            positions[index].x = vertex(0);
            positions[index].y = vertex(1);
            positions[index].z = vertex(2);
            it->next();
        }
    }
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = this->domainDimension();
    const int worldDim = this->grid()->dimWorld();
    positions.resize(m_flatLocalDofCount);

    const IndexSet& indexSet = m_view->indexSet();
    int elementCount = m_view->entityCount(0);

    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    arma::Mat<CoordinateType> corners;
    // Here we depend on the iterator returning elements in order of
    // increasing index
    size_t flatLdofIndex = 0;
    while (!it->finished())
    {
        const Entity<0>& e = it->entity();
        e.geometry().getCorners(corners);
        for (size_t v = 0; v < corners.n_cols; ++v) {
            positions[flatLdofIndex].x = corners(0, v);
            positions[flatLdofIndex].y = corners(1, v);
            positions[flatLdofIndex].z = gridDim == 2 ? corners(2, v) : 0.;
            ++flatLdofIndex;
        }
        it->next();
    }
    assert(flatLdofIndex == m_flatLocalDofCount);
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getGlobalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
   const int gridDim = this->domainDimension();
   const int globalDofCount_ = globalDofCount();
   bboxes.resize(globalDofCount_);

   const IndexSet& indexSet = m_view->indexSet();
   BoundingBox<CoordinateType> model;
   model.lbound.x = std::numeric_limits<CoordinateType>::max();
   model.lbound.y = std::numeric_limits<CoordinateType>::max();
   model.lbound.z = std::numeric_limits<CoordinateType>::max();
   model.ubound.x = -std::numeric_limits<CoordinateType>::max();
   model.ubound.y = -std::numeric_limits<CoordinateType>::max();
   model.ubound.z = -std::numeric_limits<CoordinateType>::max();
   std::fill(bboxes.begin(), bboxes.end(), model);

   std::vector<int> vertexIndices;
   arma::Mat<CoordinateType> corners;

   if (gridDim != 2)
       throw std::runtime_error("PiecewiseLinearContinuousScalarSpace::"
                                "getGlobalDofBoundingBoxes(): so far "
                                "implemented only for two-dimensional grids");

   std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
   while (!it->finished())
   {
       const Entity<0>& e = it->entity();
       const Geometry& geo = e.geometry();

       geo.getCorners(corners);
       const size_t cornerCount = corners.n_cols;
       vertexIndices.resize(cornerCount);
       for (size_t i = 0; i < cornerCount; ++i) {
           int index = indexSet.subEntityIndex(e, i, gridDim);
           bboxes[index].reference.x = corners(0, i);
           bboxes[index].reference.y = corners(1, i);
           bboxes[index].reference.z = corners(2, i);
           vertexIndices[i] = index;
       }
       for (size_t i = 0; i < cornerCount; ++i)
           for (size_t j = 0; j < cornerCount; ++j) {
               bboxes[vertexIndices[i]].lbound.x =
                   std::min(bboxes[vertexIndices[i]].lbound.x, corners(0, j));
               bboxes[vertexIndices[i]].lbound.y =
                   std::min(bboxes[vertexIndices[i]].lbound.y, corners(1, j));
               bboxes[vertexIndices[i]].lbound.z =
                   std::min(bboxes[vertexIndices[i]].lbound.z, corners(2, j));
               bboxes[vertexIndices[i]].ubound.x =
                   std::max(bboxes[vertexIndices[i]].ubound.x, corners(0, j));
               bboxes[vertexIndices[i]].ubound.y =
                   std::max(bboxes[vertexIndices[i]].ubound.y, corners(1, j));
               bboxes[vertexIndices[i]].ubound.z =
                   std::max(bboxes[vertexIndices[i]].ubound.z, corners(2, j));
           }
       it->next();
   }

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
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
getFlatLocalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
   const int gridDim = this->domainDimension();
   bboxes.resize(m_flatLocalDofCount);

   const IndexSet& indexSet = m_view->indexSet();
   BoundingBox<CoordinateType> model;
   model.lbound.x = std::numeric_limits<CoordinateType>::max();
   model.lbound.y = std::numeric_limits<CoordinateType>::max();
   model.lbound.z = std::numeric_limits<CoordinateType>::max();
   model.ubound.x = -std::numeric_limits<CoordinateType>::max();
   model.ubound.y = -std::numeric_limits<CoordinateType>::max();
   model.ubound.z = -std::numeric_limits<CoordinateType>::max();
   std::fill(bboxes.begin(), bboxes.end(), model);

   std::vector<int> vertexIndices;
   arma::Mat<CoordinateType> corners;

   if (gridDim != 2)
       throw std::runtime_error("PiecewiseLinearContinuousScalarSpace::"
                                "getFlatLocalDofBoundingBoxes(): so far "
                                "implemented only for two-dimensional grids");

   std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
   int flatLocalDof = 0;
   while (!it->finished())
   {
       const Entity<0>& e = it->entity();
       const Geometry& geo = e.geometry();

       geo.getCorners(corners);
       const size_t cornerCount = corners.n_cols;
       vertexIndices.resize(cornerCount);
       for (size_t i = 0; i < cornerCount; ++i) {
           bboxes[flatLocalDof].reference.x = corners(0, i);
           bboxes[flatLocalDof].reference.y = corners(1, i);
           bboxes[flatLocalDof].reference.z = corners(2, i);
           for (size_t j = 0; j < cornerCount; ++j) {
               bboxes[flatLocalDof].lbound.x =
                   std::min(bboxes[flatLocalDof].lbound.x, corners(0, j));
               bboxes[flatLocalDof].lbound.y =
                   std::min(bboxes[flatLocalDof].lbound.y, corners(1, j));
               bboxes[flatLocalDof].lbound.z =
                   std::min(bboxes[flatLocalDof].lbound.z, corners(2, j));
               bboxes[flatLocalDof].ubound.x =
                   std::max(bboxes[flatLocalDof].ubound.x, corners(0, j));
               bboxes[flatLocalDof].ubound.y =
                   std::max(bboxes[flatLocalDof].ubound.y, corners(1, j));
               bboxes[flatLocalDof].ubound.z =
                   std::max(bboxes[flatLocalDof].ubound.z, corners(2, j));
           }
           ++flatLocalDof;
       }
       it->next();
   }
   assert(flatLocalDof == m_flatLocalDofCount);

#ifndef NDEBUG
   std::vector<Point3D<CoordinateType> > positions;
   getFlatLocalDofPositions(positions);
   for (size_t i = 0; i < m_flatLocalDofCount; ++i) {
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
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getGlobalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    const int gridDim = this->domainDimension();
    const int globalDofCount_ = globalDofCount();
    const int worldDim = this->grid()->dimWorld();
    normals.resize(globalDofCount_);

    const IndexSet& indexSet = m_view->indexSet();
    int elementCount = m_view->entityCount(0);

    arma::Mat<CoordinateType> elementNormals(worldDim, elementCount);
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    arma::Col<CoordinateType> center(gridDim);
    center.fill(0.5);
    arma::Col<CoordinateType> normal;
    while (!it->finished())
    {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        e.geometry().getNormals(center, normal);

        for (int dim = 0; dim < worldDim; ++dim)
            elementNormals(dim, index) = normal(dim);
        it->next();
    }

    if (gridDim == 1)
        for (size_t g = 0; g < globalDofCount_; ++g) {
            normals[g].x = 0.;
            normals[g].y = 0.;
            for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
                normals[g].x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
                normals[g].y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
            }
            normals[g].x /= m_global2localDofs[g].size();
            normals[g].y /= m_global2localDofs[g].size();
        }
    else // gridDim == 2
        for (size_t g = 0; g < globalDofCount_; ++g) {
            normals[g].x = 0.;
            normals[g].y = 0.;
            for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
                normals[g].x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
                normals[g].y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
                normals[g].z += elementNormals(2, m_global2localDofs[g][l].entityIndex);
            }
            normals[g].x /= m_global2localDofs[g].size();
            normals[g].y /= m_global2localDofs[g].size();
            normals[g].z /= m_global2localDofs[g].size();
        }
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    const int gridDim = this->domainDimension();
    const int worldDim = this->grid()->dimWorld();
    normals.resize(m_flatLocalDofCount);

    const IndexSet& indexSet = m_view->indexSet();
    int elementCount = m_view->entityCount(0);

    arma::Mat<CoordinateType> elementNormals(worldDim, elementCount);
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    arma::Col<CoordinateType> center(gridDim);
    center.fill(0.5);
    arma::Col<CoordinateType> normal;
    while (!it->finished())
    {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        e.geometry().getNormals(center, normal);

        for (int dim = 0; dim < worldDim; ++dim)
            elementNormals(dim, index) = normal(dim);
        it->next();
    }

    size_t flatLdofIndex = 0;
    if (gridDim == 1)
        for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
            for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
                normals[flatLdofIndex].x = elementNormals(0, e);
                normals[flatLdofIndex].y = elementNormals(1, e);
                normals[flatLdofIndex].z = 0.;
                ++flatLdofIndex;
            }
        }
    else // gridDim == 2
        for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
            for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
                normals[flatLdofIndex].x = elementNormals(0, e);
                normals[flatLdofIndex].y = elementNormals(1, e);
                normals[flatLdofIndex].z = elementNormals(2, e);
                ++flatLdofIndex;
            }
        }
    assert(flatLdofIndex == m_flatLocalDofCount);
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs) const
{
    dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::dumpClusterIdsEx(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs,
        DofType dofType) const
{
    if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "dumpClusterIds(): invalid DOF type");
    const size_t idCount = clusterIdsOfDofs.size();
    if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
            (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "dumpClusterIds(): incorrect dimension");

    std::auto_ptr<GridView> view = this->grid()->leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();
    if (dofType == GLOBAL_DOFS) {
        arma::Row<double> data(idCount);
        for (size_t i = 0; i < idCount; ++i)
            data(i) = clusterIdsOfDofs[i];
        vtkWriter->addVertexData(data, "ids");
        vtkWriter->write(fileName);
    } else {
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
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseLinearContinuousScalarSpace);

} // namespace Bempp
