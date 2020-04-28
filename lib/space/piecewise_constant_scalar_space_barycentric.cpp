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

#include "piecewise_constant_scalar_space_barycentric.hpp"
#include "piecewise_constant_scalar_space.hpp"

#include "adaptive_space.hpp"
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

namespace Bempp {

namespace {

template <typename BasisFunctionType>
class ConstantBarycentricSpaceFactory
    : public SpaceFactory<BasisFunctionType> {
public:
  shared_ptr<Space<BasisFunctionType>>
  create(const shared_ptr<const Grid> &grid,
         const GridSegment &segment) const override {

    return shared_ptr<Space<BasisFunctionType>>(
        new PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>(
            grid, segment));
  }
};
}


template <typename BasisFunctionType>
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    PiecewiseConstantScalarSpaceBarycentric(
        const shared_ptr<const Grid> &grid)
    : ScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(GridSegment::wholeGrid(*grid)), m_strictlyOnSegment(false),
      m_putDofsOnBoundaries(true),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    PiecewiseConstantScalarSpaceBarycentric(
        const shared_ptr<const Grid> &grid, const GridSegment &segment,
        bool strictlyOnSegment, bool putDofsOnBoundaries)
    : ScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(segment), m_strictlyOnSegment(strictlyOnSegment),
      m_putDofsOnBoundaries(putDofsOnBoundaries),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::initialize() {
  const int gridDim = this->grid()->dim();
  if (gridDim != 1 && gridDim != 2)
    throw std::invalid_argument(
        "PiecewiseConstantScalarSpaceBarycentric::initialize(): "
        "only 1- and 2-dimensional grids are supported");
  m_view = this->grid()->leafView();
  assignDofsImpl();
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {
  if (!m_discontinuousSpace) {
    tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
    typedef PiecewiseConstantScalarSpace<BasisFunctionType> DiscontinuousSpace;
    if (!m_discontinuousSpace)
      m_discontinuousSpace.reset(
          new DiscontinuousSpace(m_originalGrid->barycentricGrid()));
  }
  return m_discontinuousSpace;
}

template <typename BasisFunctionType>
bool PiecewiseConstantScalarSpaceBarycentric<
    BasisFunctionType>::isDiscontinuous() const {
  return false;
}

template <typename BasisFunctionType>
int PiecewiseConstantScalarSpaceBarycentric<
    BasisFunctionType>::domainDimension() const {
  return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseConstantScalarSpaceBarycentric<
    BasisFunctionType>::codomainDimension() const {
  return 1;
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  return m_shapeset;
}

template <typename BasisFunctionType>
ElementVariant
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::elementVariant(
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
bool PiecewiseConstantScalarSpaceBarycentric<
    BasisFunctionType>::spaceIsCompatible(const Space<BasisFunctionType> &other)
    const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else {
    if (other.spaceIdentifier() == PIECEWISE_CONSTANT_SCALAR) {
      // Check if this grid is a barycentric representation of the other grid
      return this->grid()->isBarycentricRepresentationOf(*other.grid());
    } else {
      return false;
    }
  }
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::barycentricSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {

  if (self.get() != this)
    throw std::invalid_argument(
        "PiecewiseConstantScalarSpaceBarycentric::barycentricSpace(): "
        "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<
    BasisFunctionType>::setElementVariant(const Entity<0> &element,
                                          ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("PiecewiseConstantScalarSpaceBarycentric::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::assignDofsImpl() {
  const int gridDim = this->domainDimension();

  std::unique_ptr<GridView> coarseView = m_originalGrid->leafView();
  const IndexSet &index = coarseView->indexSet();
  int elementCount = this->gridView().entityCount(0);
  int faceCountCoarseGrid = coarseView->entityCount(0);

  std::vector<int> globalDofIndices(faceCountCoarseGrid, 0);
  int globalDofCount_ = 0;
  for (int faceIndex = 0; faceIndex < faceCountCoarseGrid; ++faceIndex)
    if(m_segment.contains(0,faceIndex))
      acc(globalDofIndices, faceIndex) = globalDofCount_++;
    else
      acc(globalDofIndices, faceIndex) = -1;
  int flatLocalDofCount_ = 0;

  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_global2localDofs.clear();
  m_global2localDofs.resize(globalDofCount_);

  for (std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();
       !it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    int ent0Number = index.subEntityIndex(entity, 0, 0);
    const EntityIndex faceIndex = index.subEntityIndex(entity, 0, 0);
    for (int i = 0; i != 6; ++i) {
      EntityIndex sonIndex = m_sonMap(ent0Number, i);

      std::vector<GlobalDofIndex> &globalDof =
          acc(m_local2globalDofs, sonIndex);
      globalDof.resize(1);

      GlobalDofIndex globalDofIndex;
      globalDofIndex = acc(globalDofIndices, faceIndex);
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
size_t
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::globalDofCount()
    const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t
PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::flatLocalDofCount()
    const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = this->gridView().elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                     std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i)
    localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofInterpolationPoints(Matrix<CoordinateType> &points) const {
  SpaceHelper<BasisFunctionType>::
      getGlobalDofInterpolationPoints_defaultImplementation(*this, points);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getNormalsAtGlobalDofInterpolationPoints(
        Matrix<CoordinateType> &normals) const {
  SpaceHelper<BasisFunctionType>::
      getNormalsAtGlobalDofInterpolationPoints_defaultImplementation(*this,
                                                                     normals);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getGlobalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  SpaceHelper<BasisFunctionType>::
      getGlobalDofBoundingBoxes_defaultImplementation(
          this->gridView(), m_global2localDofs, bboxes);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  getGlobalDofBoundingBoxes(bboxes);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  SpaceHelper<BasisFunctionType>::getGlobalDofNormals_defaultImplementation(
      this->gridView(), m_global2localDofs, normals);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
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
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::dumpClusterIds(
    const char *fileName,
    const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpaceBarycentric<BasisFunctionType>::
    dumpClusterIdsEx(const char *fileName,
                     const std::vector<unsigned int> &clusterIdsOfGlobalDofs,
                     DofType dofType) const {
  if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
    throw std::invalid_argument("PiecewiseConstantScalarSpaceBarycentric::"
                                "dumpClusterIds(): invalid DOF type");
  const size_t idCount = clusterIdsOfGlobalDofs.size();
  if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
      (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
    throw std::invalid_argument(
        "PiecewiseConstantScalarSpaceBarycentric::dumpClusterIds(): "
        "clusterIds has incorrect length");

  std::unique_ptr<GridView> view = this->grid()->leafView();
  std::unique_ptr<VtkWriter> vtkWriter = view->vtkWriter();
  Matrix<double> data(1, idCount);
  for (size_t i = 0; i < idCount; ++i)
    data(0, i) = clusterIdsOfGlobalDofs[i];
  vtkWriter->addCellData(data, "ids");
  vtkWriter->write(fileName);
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpaceBarycentric(
    const shared_ptr<const Grid> &grid) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new ConstantBarycentricSpaceFactory<BasisFunctionType>());
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpaceBarycentric(
    const shared_ptr<const Grid> &grid, const std::vector<int> &domains) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new ConstantBarycentricSpaceFactory<BasisFunctionType>());
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, domains, false));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseConstantScalarSpaceBarycentric(
    const shared_ptr<const Grid> &grid, int domain) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new ConstantBarycentricSpaceFactory<BasisFunctionType>());
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, std::vector<int>({domain}), false));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)                                      \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseConstantScalarSpaceBarycentric<BASIS>(              \
      const shared_ptr<const Grid> &);                                         \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseConstantScalarSpaceBarycentric<BASIS>(              \
      const shared_ptr<const Grid> &, const std::vector<int> &);   \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseConstantScalarSpaceBarycentric<BASIS>(              \
      const shared_ptr<const Grid> &, int)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewiseConstantScalarSpaceBarycentric);

} // namespace Bempp
