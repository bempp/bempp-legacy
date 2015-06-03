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

#include "piecewise_constant_dual_grid_discontinuous_scalar_space.hpp"
#include "piecewise_constant_scalar_space.hpp"

#include "space_helper.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/acc.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../common/eigen_support.hpp"

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

namespace Bempp {

template <typename BasisFunctionType>
PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    PiecewiseConstantDualGridDiscontinuousScalarSpace(
        const shared_ptr<const Grid> &grid)
    : ScalarSpace<BasisFunctionType>(grid->barycentricGrid()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::~PiecewiseConstantDualGridDiscontinuousScalarSpace() {}

template <typename BasisFunctionType>
int PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::domainDimension() const {
  return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::codomainDimension() const {
  return 1;
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  return m_basis;
}

template <typename BasisFunctionType>
ElementVariant PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::elementVariant(const Entity<0> &element) const {
  GeometryType type = element.type();
  if (type.dim() == 1)
    return 2;
  if (type.isTriangle())
    return 3;
  else
    return 4;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::setElementVariant(const Entity<0> &element,
                                          ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("PiecewiseConstantScalarSpace::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::initialize() {
  const int gridDim = this->grid()->dim();
  if (gridDim != 1 && gridDim != 2)
    throw std::invalid_argument(
        "PiecewiseConstantDualGridDiscontinuousScalarSpace::initialize(): "
        "only 1- and 2-dimensional grids are supported");
  assignDofsImpl();
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {
  if (self.get() != this)
    throw std::invalid_argument("PiecewiseConstantDualGridDiscontinuousScalarSp"
                                "ace::discontinuousSpace(): "
                                "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
bool PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::spaceIsCompatible(const Space<BasisFunctionType> &other)
    const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else
    return false;
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    barycentricSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {

  if (self.get() != this)
    throw std::invalid_argument("PiecewiseConstantDualGridDiscontinuousScalarSp"
                                "ace::barycentricSpace(): "
                                "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
bool PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::isDiscontinuous() const {
  return true;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::assignDofsImpl() {

  const GridView &view = this->gridView();
  const Mapper &elementMapper = view.elementMapper();
  size_t elementCount = view.entityCount(0);

  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_global2localDofs.clear();
  m_global2localDofs.resize(elementCount);

  // Iterate over elements of the coarse grid
  std::unique_ptr<EntityIterator<0>> it = view.entityIterator<0>();
  int flatLocalDofCount_ = 0;
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    EntityIndex elementIndex = elementMapper.entityIndex(element);
    acc(m_local2globalDofs, elementIndex).push_back(elementIndex);
    acc(m_global2localDofs, elementIndex).push_back(LocalDof(elementIndex, 0));
    ++flatLocalDofCount_;
    it->next();
  }

  // Initialize the container mapping the flat local dof indices to
  // local dof indices
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      flatLocalDofCount_, m_local2globalDofs, m_flatLocal2localDofs);
}

template <typename BasisFunctionType>
size_t PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::globalDofCount() const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::flatLocalDofCount() const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::getGlobalDofs(const Entity<0> &element,
                                      std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = this->gridView().elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                     std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i)
    localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getGlobalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  SpaceHelper<BasisFunctionType>::
      getGlobalDofBoundingBoxes_defaultImplementation(
          this->gridView(), m_global2localDofs, bboxes);
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  // TODO: extract this loop into a private function
  const IndexSet &indexSet = this->gridView().indexSet();
  const int elementCount = this->gridView().entityCount(0);

  std::vector<Matrix<CoordinateType>> elementCorners(elementCount);
  std::unique_ptr<EntityIterator<0>> it =
      this->gridView().template entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    const Geometry &geo = e.geometry();
    geo.getCorners(acc(elementCorners, index));
    it->next();
  }

  BoundingBox<CoordinateType> model;
  const CoordinateType maxCoord = std::numeric_limits<CoordinateType>::max();
  model.lbound.x = model.lbound.y = model.lbound.z = maxCoord;
  model.ubound.x = model.ubound.y = model.ubound.z = -maxCoord;

  const int flatLocalDofCount_ = m_flatLocal2localDofs.size();
  bboxes.resize(flatLocalDofCount_, model);
  for (int i = 0; i < flatLocalDofCount_; ++i) {
    const LocalDof &localDof = acc(m_flatLocal2localDofs, i);
    BoundingBox<CoordinateType> &bbox = acc(bboxes, i);
    extendBoundingBox(bbox, acc(elementCorners, localDof.entityIndex));
    setBoundingBoxReference<CoordinateType>(
        bbox, acc(elementCorners, localDof.entityIndex).col(localDof.dofIndex));
  }

#ifndef NDEBUG
  const int globalDofCount_ = globalDofCount();
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
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  SpaceHelper<BasisFunctionType>::getGlobalDofNormals_defaultImplementation(
      this->gridView(), m_global2localDofs, normals);
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = this->domainDimension();
  const int worldDim = this->grid()->dimWorld();
  normals.resize(m_flatLocal2localDofs.size());

  const IndexSet &indexSet = this->gridView().indexSet();
  int elementCount = this->gridView().entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it =
      this->gridView().template entityIterator<0>();
  Vector<CoordinateType> center(gridDim);
  center.fill(0.5);
  Matrix<CoordinateType> normal;
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getNormals(center, normal);

    for (int dim = 0; dim < worldDim; ++dim)
      elementNormals(dim, index) = normal(dim, 0);
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
void PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::dumpClusterIds(const char *fileName,
                                       const std::vector<unsigned int> &
                                           clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridDiscontinuousScalarSpace<
    BasisFunctionType>::dumpClusterIdsEx(const char *fileName,
                                         const std::vector<unsigned int> &
                                             clusterIdsOfDofs,
                                         DofType dofType) const {
  // Note: this will probably only work for spaces on full grids
  // (not on segments)

  if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
    throw std::invalid_argument(
        "PiecewiseConstantDualGridDiscontinuousScalarSpaceScalarSpace::"
        "dumpClusterIds(): invalid DOF type");
  const size_t idCount = clusterIdsOfDofs.size();
  if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
      (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
    throw std::invalid_argument(
        "PiecewiseConstantDualGridDiscontinuousScalarSpaceScalarSpace::"
        "dumpClusterIds(): incorrect dimension");

  std::unique_ptr<VtkWriter> vtkWriter = this->gridView().vtkWriter();
  if (dofType == GLOBAL_DOFS) {
    Matrix<double> data(1, idCount);
    for (size_t i = 0; i < idCount; ++i)
      data(0, i) = clusterIdsOfDofs[i];
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
  } else {
    Matrix<double> data(idCount, globalDofCount());
    data.setZero();
    size_t row = 0;
    for (size_t id = 0; id < idCount; ++id) {
      bool exists = false;
      for (size_t fldof = 0; fldof < idCount; ++fldof) {
        if (clusterIdsOfDofs[fldof] == id) {
          LocalDof ldof = m_flatLocal2localDofs[fldof];
          GlobalDofIndex gdof =
              m_local2globalDofs[ldof.entityIndex][ldof.dofIndex];
          data(row, gdof) = 1;
          exists = true;
        }
      }
      if (!exists)
        // Remove the row
        eigenRemoveRowFromMatrix(data, row);

      else
        ++row;
    }
    std::cout << "about to write" << std::endl;
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
  }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewiseConstantDualGridDiscontinuousScalarSpace);

} // namespace Bempp
