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

#include "piecewise_linear_dual_grid_continuous_scalar_space.hpp"

#include "piecewise_linear_continuous_scalar_space.hpp"
#include "space_helper.hpp"

#include "../assembly/discrete_boundary_operator.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
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
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

namespace Bempp {

namespace {

template <typename BasisFunctionType>
class LinearDualGridContinuousSpaceFactory
    : public SpaceFactory<BasisFunctionType> {
public:
  LinearDualGridContinuousSpaceFactory(bool strictlyOnSegment, bool closed)
      : m_strictlyOnSegment(strictlyOnSegment), m_closed(closed) {}
  shared_ptr<Space<BasisFunctionType>>
  create(const shared_ptr<const Grid> &grid,
         const GridSegment &segment) const override {

    return shared_ptr<Space<BasisFunctionType>>(
        new PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>(
            grid, segment, m_strictlyOnSegment, m_closed));
  }
private:
  bool m_strictlyOnSegment;
  bool m_closed;
};
}

template <typename BasisFunctionType>
PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    PiecewiseLinearDualGridContinuousScalarSpace(
        const shared_ptr<const Grid> &grid)
    : PiecewiseLinearScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(GridSegment::wholeGrid(*grid)), m_strictlyOnSegment(false),
      m_putDofsOnBoundaries(true),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    PiecewiseLinearDualGridContinuousScalarSpace(
        const shared_ptr<const Grid> &grid, const GridSegment &segment,
        bool strictlyOnSegment, bool putDofsOnBoundaries)
    : PiecewiseLinearScalarSpace<BasisFunctionType>(grid->barycentricGrid()),
      m_segment(segment), m_strictlyOnSegment(strictlyOnSegment),
      m_putDofsOnBoundaries(putDofsOnBoundaries),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::~PiecewiseLinearDualGridContinuousScalarSpace() {}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::initialize() {
  const int gridDim = this->grid()->dim();
  if (gridDim != 1 && gridDim != 2)
    throw std::invalid_argument(
        "PiecewiseLinearDualGridContinuousScalarSpace::initialize(): "
        "only 1- and 2-dimensional grids are supported");
  m_view = this->grid()->leafView();
  assignDofsImpl();
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {
    throw std::runtime_error("PiecewiseLinearDualGridContinuousScalarSpace::discontinuousSpace(): "
                             "Not yet implemented.");
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    barycentricSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {

  if (self.get() != this)
    throw std::invalid_argument(
        "PiecewiseLinearDualGridContinuousScalarSpace::barycentricSpace(): "
        "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
int PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::domainDimension() const {
  return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::codomainDimension() const {
  return 1;
}

template <typename BasisFunctionType>
bool PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::isDiscontinuous() const {
  return false;
}

template <typename BasisFunctionType>
bool PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::spaceIsCompatible(const Space<BasisFunctionType> &other)
    const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else {
    if (other.spaceIdentifier() == PIECEWISE_LINEAR_CONTINUOUS_SCALAR) {
      // Check if this grid is a barycentric representation of the other grid
      return this->grid()->isBarycentricRepresentationOf(*other.grid());
    } else {
      return false;
    }
  }
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::assignDofsImpl() {
  // Set up useful numbers, maps, etc.
  std::unique_ptr<GridView> coarseView = m_originalGrid->leafView();
  const int faceCountCoarseGrid = coarseView->entityCount(0);
  const int edgeCountCoarseGrid = coarseView->entityCount(1);
  const int vertexCountCoarseGrid = coarseView->entityCount(2);
  const int faceCountFineGrid = m_view->entityCount(0);
  const int edgeCountFineGrid = m_view->entityCount(1);
  const int vertexCountFineGrid = m_view->entityCount(2);
  const IndexSet &index = coarseView->indexSet();
  const IndexSet &bindex = m_view->indexSet();

  // Set up useful numbers, maps, etc.
  std::vector<bool> fineFaceInSegment(faceCountFineGrid, false);
  std::vector<int> fineFaceCountAroundCoarseVertex(vertexCountCoarseGrid, 0);

  for(std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();!it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = index.subEntityIndex(entity, 0, 0);
    for(int i=0;i!=3;++i)
      fineFaceCountAroundCoarseVertex[index.subEntityIndex(entity, i, 2)] += 2;
    for(int i=0;i!=6;++i)
      fineFaceInSegment[m_sonMap(ent0Number, i)] = m_segment.contains(0, ent0Number);
  }

  std::vector<bool> fullyInSegment;
  fullyInSegment.resize(faceCountCoarseGrid, true);
  std::vector<bool> vertexFullyInSegment;
  vertexFullyInSegment.resize(vertexCountCoarseGrid, true);
  std::vector<bool> isAttachedToSegment;
  isAttachedToSegment.resize(faceCountCoarseGrid, false);
  std::vector<bool> vertexAttachedToSegment;
  vertexAttachedToSegment.resize(vertexCountCoarseGrid, false);

  Matrix<int> verticesOfFineFace;
  verticesOfFineFace.resize(faceCountFineGrid,3);
  for(std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();!it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = bindex.entityIndex(entity);
    for(int i=0;i!=3;++i)
      verticesOfFineFace(ent0Number,i) = bindex.subEntityIndex(entity, i, 2);
  }

  for(std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();!it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = index.subEntityIndex(entity, 0, 0);
    for(int i=0; i!=3; ++i) {
      const int ent2Number = index.subEntityIndex(entity, i, 2);
      if(m_segment.contains(0,ent0Number))
        vertexAttachedToSegment[ent2Number] = true;
      else
        vertexFullyInSegment[ent2Number] = false;
    }
  }
  for(std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();!it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = index.entityIndex(entity);
    bool att = false;
    for(int i=0; i!=3; ++i) {
      const int ent2Number = index.subEntityIndex(entity, i, 2);
      if(!vertexFullyInSegment[ent2Number])
        fullyInSegment[ent0Number] = false;
      if(vertexAttachedToSegment[ent2Number])
        isAttachedToSegment[ent0Number] = true;
    }
  }

  // Assign Dofs to Faces
  std::vector<int> globalDofsOfFaces;
  globalDofsOfFaces.resize(faceCountCoarseGrid);
  int globalDofCount_ = 0;
  for (int i = 0; i != faceCountCoarseGrid; ++i) {
    int &globalDofOfFace = acc(globalDofsOfFaces, i);

    if(fullyInSegment[i])
      globalDofOfFace = globalDofCount_++;
    else if (m_putDofsOnBoundaries && isAttachedToSegment[i])
      globalDofOfFace = globalDofCount_++;
    else
      globalDofOfFace = -1;
  }


  // Map out dofs
  std::vector<std::vector<int>> dofsAtFineVertex;
  dofsAtFineVertex.resize(vertexCountFineGrid);
  std::vector<std::vector<BasisFunctionType>> dofsWeightsAtFineVertex;
  dofsWeightsAtFineVertex.resize(vertexCountFineGrid);

  for(std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();!it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = index.entityIndex(entity);
    const int son0 = m_sonMap(ent0Number, 0);
    const int son2 = m_sonMap(ent0Number, 2);
    const int son4 = m_sonMap(ent0Number, 4);
    const int dof = globalDofsOfFaces[ent0Number];
    if(dof != -1){
      dofsAtFineVertex[verticesOfFineFace(son0,1)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son0,1)].push_back(1.0);
      dofsAtFineVertex[verticesOfFineFace(son0,2)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son0,2)].push_back(0.5);
      dofsAtFineVertex[verticesOfFineFace(son2,2)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son2,2)].push_back(0.5);
      dofsAtFineVertex[verticesOfFineFace(son4,2)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son4,2)].push_back(0.5);
      dofsAtFineVertex[verticesOfFineFace(son0,0)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son0,0)].push_back(1./fineFaceCountAroundCoarseVertex[index.subEntityIndex(entity,0,2)]);
      dofsAtFineVertex[verticesOfFineFace(son2,0)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son2,0)].push_back(1./fineFaceCountAroundCoarseVertex[index.subEntityIndex(entity,1,2)]);
      dofsAtFineVertex[verticesOfFineFace(son4,0)].push_back(dof);
      dofsWeightsAtFineVertex[verticesOfFineFace(son4,0)].push_back(1./fineFaceCountAroundCoarseVertex[index.subEntityIndex(entity,2,2)]);
    }
  }





  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(faceCountFineGrid);
  m_global2localDofs.clear();
  m_global2localDofs.resize(globalDofCount_);
  m_fineFaceCoeffs.clear();
  m_fineFaceCoeffs.resize(faceCountFineGrid);
  size_t flatLocalDofCount_ = 0;

  // Initialise bounding-box caches
  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  m_globalDofBoundingBoxes.resize(globalDofCount_, model);
  m_elementShapesets.resize(faceCountFineGrid);

  // Set up coefficients for shapesets
  for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
       !it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = bindex.entityIndex(entity);
    if(!m_strictlyOnSegment || fineFaceInSegment[ent0Number]){
      for(int v=0;v!=3;++v){
        const int ent2Number = bindex.subEntityIndex(entity, v, 2);
        Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[ent0Number];
        for(int i=0;i!=dofsAtFineVertex[ent2Number].size();++i){
          const int glDof = dofsAtFineVertex[ent2Number][i];
          bool done = false;
          for(int d=0;d<acc(m_local2globalDofs, ent0Number).size() && !done;++d){
            if(acc(m_local2globalDofs, ent0Number)[d] == glDof){
              ffCoeff(v, d) += dofsWeightsAtFineVertex[ent2Number][i];
              done = true;
            }
          }
          if(!done){
            const int M = ffCoeff.cols();
            ffCoeff.conservativeResize(3, M+1);
            for(int j=0;j!=3;++j)
              ffCoeff(j, M) = 0.;
            ffCoeff(v, M) = dofsWeightsAtFineVertex[ent2Number][i];

            acc(m_local2globalDofs, ent0Number).push_back(glDof);
            if(glDof!=-1){
              acc(m_global2localDofs, glDof).push_back(LocalDof(ent0Number, M+i));
              ++flatLocalDofCount_;
            }
          }
        }
      }
    }
  }


  // Bounding boxes
  for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
       !it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    const int ent0Number = bindex.entityIndex(entity);
    Matrix<CoordinateType> vertices;
    const Geometry &geo = entity.geometry();
    geo.getCorners(vertices);

    for (int j = 0; j != m_local2globalDofs[ent0Number].size(); ++j) {
      int glDof = m_local2globalDofs[ent0Number][j];
      extendBoundingBox(acc(m_globalDofBoundingBoxes, glDof), vertices);
      setBoundingBoxReference<CoordinateType>(
          acc(m_globalDofBoundingBoxes, glDof),
          0.5 * (vertices.col(0) + vertices.col(1)));
    }
  }

  for (std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
       !it->finished(); it->next()) {
    const Entity<0> &entity = it->entity();
    int ent0Number = bindex.entityIndex(entity);
    Matrix<BasisFunctionType> &ffCoeff = m_fineFaceCoeffs[ent0Number];
    if (ffCoeff.cols() == 0) {
      m_local2globalDofs[ent0Number].push_back(-1);
      ffCoeff.conservativeResize(3, 1);
      ffCoeff(0, 0) = 0;
      ffCoeff(1, 0) = 0;
      ffCoeff(2, 0) = 0;
    }
    m_elementShapesets[ent0Number] = Shapeset(ffCoeff);
  }
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      flatLocalDofCount_, m_local2globalDofs, m_flatLocal2localDofs);

}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  const Mapper &elementMapper = m_view->elementMapper();
  int index = elementMapper.entityIndex(element);
  return m_elementShapesets[index];
}

template <typename BasisFunctionType>
ElementVariant PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::elementVariant(const Entity<0> &element) const {
  GeometryType type = element.type();
  if (type.isLine())
    return 2;
  else if (type.isTriangle())
    return 3;
  else if (type.isQuadrilateral())
    return 4;
  else
    throw std::runtime_error("PiecewiseLinearScalarSpace::"
                             "elementVariant(): invalid geometry type, "
                             "this shouldn't happen!");
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::setElementVariant(const Entity<0> &element,
                                          ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("PiecewiseLinearScalarSpace::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
size_t PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::globalDofCount() const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::flatLocalDofCount() const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<
    BasisFunctionType>::getGlobalDofs(const Entity<0> &element,
                                      std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = m_view->elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = acc(m_local2globalDofs, index);
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                     std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i) {
    acc(localDofs, i) = acc(m_global2localDofs, acc(globalDofs, i));
    for (size_t j = 0; j < localDofs[i].size(); ++j)
      LocalDof ldof = acc(localDofs[i], j);
  }
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  positions.resize(m_globalDofBoundingBoxes.size());
  for (size_t i = 0; i < m_globalDofBoundingBoxes.size(); ++i)
    acc(positions, i) = acc(m_globalDofBoundingBoxes, i).reference;

}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);

  positions.resize(bboxes.size());
  for (int i = 0; i < positions.size(); ++i)
    positions[i] = bboxes[i].reference;
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
    bboxes = m_globalDofBoundingBoxes;
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  const int flatLocalDofCount_ = m_flatLocal2localDofs.size();
  bboxes.resize(flatLocalDofCount_);

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  std::vector<Matrix<CoordinateType>> elementCorners(elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  while (!it->finished()) {
    const Entity<0> &e = it->entity();
    int index = indexSet.entityIndex(e);
    e.geometry().getCorners(acc(elementCorners, index));
    if (acc(elementCorners, index).cols() != 3)
      throw std::runtime_error(
          "PiecewiseLinearDualGridContinuousScalarSpace::getFlatLocalDofBoundingBoxes(): "
          "only triangular elements are supported at present");
    it->next();
  }

  size_t flatLdofIndex = 0;
  Vector<CoordinateType> dofPosition;
  for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
    for (size_t v = 0; v < acc(m_local2globalDofs, e).size(); ++v) {
      if (acc(acc(m_local2globalDofs, e), v) >= 0) { // is this LDOF used?
        const Matrix<CoordinateType> &vertices = acc(elementCorners, e);
        BoundingBox<CoordinateType> &bbox = acc(bboxes, flatLdofIndex);
        if (v == 0)
          dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
        else if (v == 1)
          dofPosition = 0.5 * (vertices.col(2) + vertices.col(0));
        else // v == 2
          dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
        extendBoundingBox(bbox, vertices);
        setBoundingBoxReference<CoordinateType>(bbox, dofPosition);
        ++flatLdofIndex;
      }
    }
  }
  assert(flatLdofIndex == flatLocalDofCount_);
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = this->domainDimension();
  const int globalDofCount_ = globalDofCount();
  const int worldDim = this->grid()->dimWorld();
  normals.resize(globalDofCount_);

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
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
    for (size_t g = 0; g < globalDofCount_; ++g) {
      normals[g].x = 0.;
      normals[g].y = 0.;
      for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
        normals[g].x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
        normals[g].y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
      }
      auto len =
          std::sqrt(normals[g].x * normals[g].x + normals[g].y * normals[g].y);
      normals[g].x /= len;
      normals[g].y /= len;
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
      auto len =
          std::sqrt(normals[g].x * normals[g].x + normals[g].y * normals[g].y +
                    normals[g].z * normals[g].z);
      normals[g].x /= len;
      normals[g].y /= len;
      normals[g].z /= len;
    }
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = this->domainDimension();
  const int worldDim = this->grid()->dimWorld();
  normals.resize(m_flatLocal2localDofs.size());

  const IndexSet &indexSet = m_view->indexSet();
  int elementCount = m_view->entityCount(0);

  Matrix<CoordinateType> elementNormals(worldDim, elementCount);
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
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
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    dumpClusterIds(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewiseLinearDualGridContinuousScalarSpace<BasisFunctionType>::
    dumpClusterIdsEx(const char *fileName,
                     const std::vector<unsigned int> &clusterIdsOfDofs,
                     DofType dofType) const {
  // Note: this will probably only work for spaces on full grids
  // (not on segments)

  if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
    throw std::invalid_argument(
        "PiecewiseLinearDualGridContinuousScalarSpace::"
        "dumpClusterIds(): invalid DOF type");
  const size_t idCount = clusterIdsOfDofs.size();
  if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
      (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
    throw std::invalid_argument(
        "PiecewiseLinearDualGridContinuousScalarSpace::"
        "dumpClusterIds(): incorrect dimension");

  std::unique_ptr<GridView> view = this->grid()->leafView();
  std::unique_ptr<VtkWriter> vtkWriter = view->vtkWriter();
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
    std::cout << "about to write" << std::endl;
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
  }
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseLinearDualGridContinuousScalarSpace(
    const shared_ptr<const Grid> &grid) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new LinearDualGridContinuousSpaceFactory<BasisFunctionType>(false, true));
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseLinearDualGridContinuousScalarSpace(
    const shared_ptr<const Grid> &grid, const std::vector<int> &domains,
    bool open, bool strictlyOnSegment) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new LinearDualGridContinuousSpaceFactory<BasisFunctionType>(strictlyOnSegment, open));
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, domains, open));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptivePiecewiseLinearDualGridContinuousScalarSpace(
    const shared_ptr<const Grid> &grid, int domain,
    bool open, bool strictlyOnSegment) {

  shared_ptr<SpaceFactory<BasisFunctionType>> factory(
      new LinearDualGridContinuousSpaceFactory<BasisFunctionType>(strictlyOnSegment, open));
  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType>(factory, grid, std::vector<int>({domain}), open));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)                                      \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseLinearDualGridContinuousScalarSpace<BASIS>(              \
      const shared_ptr<const Grid> &);                                         \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseLinearDualGridContinuousScalarSpace<BASIS>(              \
      const shared_ptr<const Grid> &, const std::vector<int> &, bool, bool);   \
  template shared_ptr<Space<BASIS>>                                            \
  adaptivePiecewiseLinearDualGridContinuousScalarSpace<BASIS>(              \
      const shared_ptr<const Grid> &, int, bool, bool)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewiseLinearDualGridContinuousScalarSpace);

} // namespace Bempp
