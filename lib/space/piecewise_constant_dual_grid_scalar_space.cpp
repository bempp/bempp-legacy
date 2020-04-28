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

#include "piecewise_constant_dual_grid_scalar_space.hpp"

#include "adaptive_space.hpp"
#include "space_helper.hpp"
#include "piecewise_constant_dual_grid_discontinuous_scalar_space.hpp"

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

namespace {

template <typename BasisFunctionType>
class PiecewiseConstantDualGridSpaceFactory
    : public SpaceFactory<BasisFunctionType> {
public:
  PiecewiseConstantDualGridSpaceFactory(bool strictlyOnSegment)
      : m_strictlyOnSegment(strictlyOnSegment) {}
  shared_ptr<Space<BasisFunctionType>>
  create(const shared_ptr<const Grid> &grid,
         const GridSegment &segment) const override {

    return shared_ptr<Space<BasisFunctionType>>(
        new PiecewiseConstantDualGridScalarSpace<BasisFunctionType>(grid,
                                                                    segment,
                                                                    m_strictlyOnSegment));
  }

private:
  bool m_strictlyOnSegment;
};
}

template <typename BasisFunctionType>
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    PiecewiseConstantDualGridScalarSpace(const shared_ptr<const Grid> &grid)
    : ScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(GridSegment::wholeGrid(*grid)), m_strictlyOnSegment(false),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}
template <typename BasisFunctionType>
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    PiecewiseConstantDualGridScalarSpace(const shared_ptr<const Grid> &grid,
                                         const GridSegment &segment,
                                         bool strictlyOnSegment)
    : ScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(segment), m_strictlyOnSegment(strictlyOnSegment),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseConstantDualGridScalarSpace<
    BasisFunctionType>::~PiecewiseConstantDualGridScalarSpace() {}

template <typename BasisFunctionType>
int PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::domainDimension()
    const {
  return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::codomainDimension()
    const {
  return 1;
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  return m_basis;
}

template <typename BasisFunctionType>
ElementVariant
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::elementVariant(
    const Entity<0> &element) const {
  GeometryType type = element.type();
  if (type.dim() == 1)
    return 2;
  if (type.isTriangle())
    return 3;
  else
    return 4;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::setElementVariant(
    const Entity<0> &element, ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("PiecewiseConstantScalarSpace::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::initialize() {
  const int gridDim = this->grid()->dim();
  if (gridDim != 1 && gridDim != 2)
    throw std::invalid_argument(
        "PiecewiseConstantDualGridScalarSpace::initialize(): "
        "only 1- and 2-dimensional grids are supported");
  assignDofsImpl();
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {
  if (!m_discontinuousSpace) {
    tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
    if (!m_discontinuousSpace)
      m_discontinuousSpace.reset(
          new PiecewiseConstantDualGridDiscontinuousScalarSpace<
              BasisFunctionType>(this->grid()));
  }
  return m_discontinuousSpace;
}

template <typename BasisFunctionType>
bool PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::spaceIsCompatible(
    const Space<BasisFunctionType> &other) const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else
    return false;
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::barycentricSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {

  if (self.get() != this)
    throw std::invalid_argument(
        "PiecewiseConstantDualGridScalarSpace::barycentricSpace(): "
        "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
bool PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::isDiscontinuous()
    const {
  return false;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::assignDofsImpl() {

  const int gridDim = this->domainDimension();

  std::unique_ptr<GridView> coarseView = m_originalGrid->leafView();
  const IndexSet &index = coarseView->indexSet();
  int elementCount = this->gridView().entityCount(0);
  int vertexCountCoarseGrid = coarseView->entityCount(gridDim);

  std::vector<int> globalDofIndices(vertexCountCoarseGrid, 0);
  m_segment.markExcludedEntities(gridDim, globalDofIndices);
  int globalDofCount_ = 0;
  for (int vertexIndex = 0; vertexIndex < vertexCountCoarseGrid; ++vertexIndex)
    if(acc(globalDofIndices, vertexIndex)==0)
      acc(globalDofIndices, vertexIndex) = globalDofCount_++;
  int flatLocalDofCount_ = 0;

  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_global2localDofs.clear();
  m_global2localDofs.resize(globalDofCount_);

  for (std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();
       !it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    int ent0Number = index.subEntityIndex(entity, 0, 0);
    bool elementContained = m_strictlyOnSegment ? m_segment.contains(0,ent0Number) : true;
    for (int i = 0; i != 6; ++i) {
      EntityIndex sonIndex = m_sonMap(ent0Number, i);

      std::vector<GlobalDofIndex> &globalDof =
          acc(m_local2globalDofs, sonIndex);
      globalDof.resize(1); // 3

      EntityIndex vertexIndex = index.subEntityIndex(entity, i / 2, gridDim);
      GlobalDofIndex globalDofIndex;
      globalDofIndex = elementContained ? acc(globalDofIndices, vertexIndex) : -1;
      //        globalDof.push_back(globalDofIndex);
      globalDof[0] = globalDofIndex;
      if (globalDofIndex >= 0) {
        acc(m_global2localDofs, globalDofIndex)
            .push_back(LocalDof(sonIndex, 0));
        ++flatLocalDofCount_;
      }
    }
  }
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      flatLocalDofCount_, m_local2globalDofs, m_flatLocal2localDofs);
}

template <typename BasisFunctionType>
size_t PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::globalDofCount()
    const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t
PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::flatLocalDofCount()
    const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = this->gridView().elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::global2localDofs(
    const std::vector<GlobalDofIndex> &globalDofs,
    std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i)
    localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getGlobalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  SpaceHelper<BasisFunctionType>::
      getGlobalDofBoundingBoxes_defaultImplementation(
          this->gridView(), m_global2localDofs, bboxes);
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
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
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  SpaceHelper<BasisFunctionType>::getGlobalDofNormals_defaultImplementation(
      this->gridView(), m_global2localDofs, normals);
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::
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
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::dumpClusterIds(
    const char *fileName,
    const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseConstantDualGridScalarSpace<BasisFunctionType>::dumpClusterIdsEx(
    const char *fileName, const std::vector<unsigned int> &clusterIdsOfDofs,
    DofType dofType) const {
  // Note: this will probably only work for spaces on full grids
  // (not on segments)

  if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
    throw std::invalid_argument(
        "PiecewiseConstantDualGridScalarSpaceScalarSpace::"
        "dumpClusterIds(): invalid DOF type");
  const size_t idCount = clusterIdsOfDofs.size();
  if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
      (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
    throw std::invalid_argument(
        "PiecewiseConstantDualGridScalarSpaceScalarSpace::"
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
    data.fill(0.);
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
        eigenRemoveRowFromMatrix(data, row);
      else
        ++row;
    }
    // std::cout << "about to write" << std::endl;
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
  }
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantDualGridScalarSpace(
    const shared_ptr<const Grid> &grid) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new PiecewiseConstantDualGridSpaceFactory<BasisFunctionType>(false));
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantDualGridScalarSpace(
    const shared_ptr<const Grid> &grid, const std::vector<int> &domains,
    bool open, bool strictlyOnSegment) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new PiecewiseConstantDualGridSpaceFactory<BasisFunctionType>(strictlyOnSegment));
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, domains, open));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantDualGridScalarSpace(
    const shared_ptr<const Grid> &grid, int domain,
    bool open, bool strictlyOnSegment) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new PiecewiseConstantDualGridSpaceFactory<BasisFunctionType>(strictlyOnSegment));
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, std::vector<int>({domain}), open));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)                                      \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseConstantDualGridScalarSpace<BASIS>(                         \
      const shared_ptr<const Grid> &);                                         \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseConstantDualGridScalarSpace<BASIS>(                         \
      const shared_ptr<const Grid> &, const std::vector<int> &, bool, bool);   \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseConstantDualGridScalarSpace<BASIS>(                         \
      const shared_ptr<const Grid> &, int, bool, bool)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewiseConstantDualGridScalarSpace);

} // namespace Bempp
