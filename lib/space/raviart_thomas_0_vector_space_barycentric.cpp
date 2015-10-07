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

#include "raviart_thomas_0_vector_space_barycentric.hpp"
#include "adaptive_space.hpp"

#include "piecewise_linear_discontinuous_scalar_space.hpp"
#include "space_helper.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/acc.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/bounding_box.hpp"
#include "../common/bounding_box_helpers.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/hdiv_function_value_functor.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
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

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::Impl {
  typedef Fiber::HdivFunctionValueFunctor<CoordinateType> TransformationFunctor;

  Impl() : transformations(TransformationFunctor()) {}

  Fiber::DefaultCollectionOfShapesetTransformations<TransformationFunctor>
      transformations;
};
/** \endcond */

template <typename BasisFunctionType>
RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::RaviartThomas0VectorSpaceBarycentric(
    const shared_ptr<const Grid> &grid, bool putDofsOnBoundaries)
    : Base(grid->barycentricGrid()), m_impl(new Impl), m_segment(GridSegment::wholeGrid(*grid)),
      m_putDofsOnBoundaries(putDofsOnBoundaries), m_dofMode(EDGE_ON_SEGMENT), m_vectorBasisType1(Shapeset::TYPE1),
      m_vectorBasisType2(Shapeset::TYPE2),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  initialize();
}

template <typename BasisFunctionType>
RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::RaviartThomas0VectorSpaceBarycentric(
    const shared_ptr<const Grid> &grid, const GridSegment &segment,
    bool putDofsOnBoundaries, int dofMode)
    : Base(grid->barycentricGrid()), m_impl(new Impl), m_segment(segment),
      m_putDofsOnBoundaries(putDofsOnBoundaries), m_dofMode(dofMode), m_vectorBasisType1(Shapeset::TYPE1),
      m_vectorBasisType2(Shapeset::TYPE2),
      m_originalGrid(grid), m_sonMap(grid->barycentricSonMap()) {
  if (!(dofMode & (EDGE_ON_SEGMENT | ELEMENT_ON_SEGMENT)))
    throw std::invalid_argument("RaviartThomas0VectorSpaceBarycentric::"
                                "RaviartThomas0VectorSpaceBarycentric(): "
                                "invalid dofMode");
  initialize();
}

template <typename BasisFunctionType>
bool RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::spaceIsCompatible(
    const Space<BasisFunctionType> &other) const {

  if (other.grid().get() == this->grid().get()) {
    return (other.spaceIdentifier() == this->spaceIdentifier());
  } else
    return false;
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::initialize() {
  if (this->grid()->dim() != 2 || this->grid()->dimWorld() != 3)
    throw std::invalid_argument("RaviartThomas0VectorSpaceBarycentric::initialize(): "
                                "grid must be 2-dimensional and embedded "
                                "in 3-dimensional space");
  m_view = this->grid()->leafView();
  std::cout << "BELLO!" << std::endl;
  assignDofsImpl();
}

template <typename BasisFunctionType>
RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::~RaviartThomas0VectorSpaceBarycentric() {}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType>> &self) const {
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
bool RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::isDiscontinuous() const {
  return false;
}

template <typename BasisFunctionType>
const typename RaviartThomas0VectorSpaceBarycentric<
    BasisFunctionType>::CollectionOfShapesetTransformations &
RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::basisFunctionValue() const {
  return m_impl->transformations;
}

template <typename BasisFunctionType>
int RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::domainDimension() const {
  return 2;
}

template <typename BasisFunctionType>
int RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::codomainDimension() const {
  return 3;
}

template <typename BasisFunctionType>
ElementVariant RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::elementVariant(
    const Entity<0> &element) const {
  GeometryType type = element.type();
  if (type.isTriangle())
    return 3;
  else if (type.isQuadrilateral())
    return 4;
  else
    throw std::runtime_error("RaviartThomas0VectorSpaceBarycentric::"
                             "elementVariant(): invalid geometry type, "
                             "this shouldn't happen!");
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::setElementVariant(
    const Entity<0> &element, ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("RaviartThomas0VectorSpaceBarycentric::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::assignDofsImpl() {
  std::cout << "1";
  int edgeCountCoarseGrid = m_originalGrid->leafView()->entityCount(1);
  int elementCount = m_view->entityCount(0);
  std::unique_ptr<GridView> coarseView = m_originalGrid->leafView();
  const IndexSet &index = coarseView->indexSet();

  std::vector<int> globalDofsOfEdges;
  globalDofsOfEdges.resize(edgeCountCoarseGrid);
  int globalDofCount_ = 0;
  for (int i = 0; i != edgeCountCoarseGrid; ++i) {
      int &globalDofOfEdge = acc(globalDofsOfEdges,i);
      globalDofOfEdge = globalDofCount_++;
  }

  std::vector<int> lowestIndicesOfElementsAdjacentToEdges(edgeCountCoarseGrid, std::numeric_limits<int>::max());

  for(std::unique_ptr<EntityIterator<0>> it = coarseView->entityIterator<0>();!it->finished();it->next()){
    for(int i=0;i!=3;++i){
      const Entity<0> &entity = it->entity();
        const int ent0Number = index.subEntityIndex(entity,0,0);

       int &lowestIndex = lowestIndicesOfElementsAdjacentToEdges[index.subEntityIndex(entity,i,1)];
       lowestIndex = std::min(ent0Number,lowestIndex);
    }
  }

  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_local2globalDofWeights.clear();
  m_local2globalDofWeights.resize(elementCount);
  m_global2localDofs.clear();
  m_global2localDofs.resize(globalDofCount_);
  size_t flatLocalDofCount = 0;

  // Initialise bounding-box caches
  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  m_globalDofBoundingBoxes.resize(globalDofCount_, model);
  m_elementIndex2Type.resize(elementCount);

  const int element2Basis[6][3] = {
                                   {0,2,1},
                                   {0,2,1},
                                   {2,0,1},
                                   {2,0,1},
                                   {1,2,0},
                                   {1,2,0}
                                  };

  int i__ = 0;
  Vector<CoordinateType> dofPosition;
  Matrix<CoordinateType> vertices;
  for(std::unique_ptr<EntityIterator<0>> it=coarseView->entityIterator<0>();!it->finished();it->next()){
    const Entity<0> &entity = it->entity();
    const Geometry &geo = entity.geometry();
    int ent0Number = index.subEntityIndex(entity,0,0);
    geo.getCorners(vertices);

    std::vector<int> edges;
    edges.resize(3);
    edges[0] = index.subEntityIndex(entity,0,1);
    edges[1] = index.subEntityIndex(entity,1,1);
    edges[2] = index.subEntityIndex(entity,2,1);
    for(int i=0;i!=6;++i){
      EntityIndex sonIndex = m_sonMap(ent0Number,i);
      if (i % 2 == 0) {
        acc(m_elementIndex2Type, sonIndex) = Shapeset::TYPE1;
      } else {
        acc(m_elementIndex2Type, sonIndex) = Shapeset::TYPE2;
      }
      std::vector<GlobalDofIndex> &globalDof = acc(m_local2globalDofs, sonIndex);
      globalDof.resize(3);

      std::vector<BasisFunctionType> &globalDofWeights = acc(m_local2globalDofWeights, sonIndex);
      globalDofWeights.resize(3);

      for(int j=0;j!=3;++j){
        const int edgeIndex = edges[element2Basis[i][j]];
        const int globalDofIndex = globalDofsOfEdges[edgeIndex];
        globalDof[j] = globalDofIndex;
        globalDofWeights[j] = acc(lowestIndicesOfElementsAdjacentToEdges, edgeIndex) == ent0Number
              ? 1.
              : -1.;
        m_global2localDofs[globalDofIndex].push_back(LocalDof(ent0Number,j));
        if(element2Basis[i][j]==0){
          dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
        } else if (element2Basis[i][j] == 1) {
          dofPosition = 0.5 * (vertices.col(2) + vertices.col(0));
        } else { // element2Basis[i][j] == 2
          dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
        }
        extendBoundingBox(acc(m_globalDofBoundingBoxes, globalDofIndex), vertices);
        setBoundingBoxReference<CoordinateType>(
            acc(m_globalDofBoundingBoxes, globalDofIndex), dofPosition);

        ++flatLocalDofCount;
      }
    }
    std::cout << i__++ << " ";
  }



/*  const Mapper &elementMapper = m_view->elementMapper();
  const IndexSet &indexSet = m_view->indexSet();

  int edgeCount = m_view->entityCount(1);
  int elementCount = m_view->entityCount(0);
  std::cout << m_view->entityCount(0) << std::endl;
  std::cout << m_originalGrid->leafView()->entityCount(0) << std::endl;
  // std::cout << "dofMode: " << m_dofMode << std::endl;

  const int edgeCodim = 1;
  const int elementCodim = 0;

  std::vector<int> lowestIndicesOfElementsAdjacentToEdges(
      edgeCount, std::numeric_limits<int>::max());
  // number of element adjacent to each edge
  std::vector<int> elementsAdjacentToEdges(edgeCount, 0);
  std::vector<bool> noAdjacentElementsAreInSegment(edgeCount, true);
  std::unique_ptr<EntityIterator<elementCodim>> it =
      m_view->entityIterator<elementCodim>();
  while (!it->finished()) {
    const Entity<elementCodim> &element = it->entity();
    const int elementIndex = indexSet.entityIndex(element);
    const bool elementContained =
        m_segment.contains(elementCodim, elementIndex);
    const int localEdgeCount = element.subEntityCount<edgeCodim>();
    for (int i = 0; i < localEdgeCount; ++i) {
      int edgeIndex = indexSet.subEntityIndex(element, i, edgeCodim);
      if (m_dofMode & EDGE_ON_SEGMENT &&
          !m_segment.contains(edgeCodim, edgeIndex))
        continue;
      ++acc(elementsAdjacentToEdges, edgeIndex);
      if (elementContained)
        acc(noAdjacentElementsAreInSegment, edgeIndex) = false;
      int &lowestIndex = acc(lowestIndicesOfElementsAdjacentToEdges, edgeIndex);
      lowestIndex = std::min(lowestIndex, elementIndex);
    }
    it->next();
  }

  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_local2globalDofWeights.clear();
  m_local2globalDofWeights.resize(elementCount);
  m_global2localDofs.clear();
  m_global2localDofs.resize(globalDofCount_);
  size_t flatLocalDofCount = 0;

  // TODO: consider calling reserve(2) for each element of m_global2localDofs

  // Initialise bounding-box caches
  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  m_globalDofBoundingBoxes.resize(globalDofCount_, model);

  // Iterate over elements
  it = m_view->entityIterator<elementCodim>();
  Matrix<CoordinateType> vertices;
  Vector<CoordinateType> dofPosition;
  while (!it->finished()) {
    const Entity<elementCodim> &element = it->entity();
    const Geometry &geo = element.geometry();
    EntityIndex elementIndex = elementMapper.entityIndex(element);
    bool elementContained = m_dofMode & ELEMENT_ON_SEGMENT
                                ? m_segment.contains(elementCodim, elementIndex)
                                : true;

    geo.getCorners(vertices);
    const int vertexCount = vertices.cols();
    const int edgeCount = vertexCount;
    if (edgeCount != 3)
      throw std::runtime_error(
          "RaviartThomas0VectorSpaceBarycentric::"
          "assignDofsImpl(): support for quadrilaterals not in place yet");

    // List of global DOF indices corresponding to the local DOFs of the
    // current element
    std::vector<GlobalDofIndex> &globalDofs =
        acc(m_local2globalDofs, elementIndex);
    // List of weights of the local DOFs residing on the current element
    std::vector<BasisFunctionType> &globalDofWeights =
        acc(m_local2globalDofWeights, elementIndex);
    globalDofs.reserve(edgeCount);
    globalDofWeights.reserve(edgeCount);
    for (int i = 0; i < edgeCount; ++i) {
      int edgeIndex = indexSet.subEntityIndex(element, i, edgeCodim);
      GlobalDofIndex globalDofIndex =
          elementContained ? acc(globalDofsOfEdges, edgeIndex) : -1;
      // Handle Dune's funny subentity indexing order. This code is
      // only valid for triangular elements.
      if (i == 0) {
        dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
      } else if (i == 1) {
        dofPosition = 0.5 * (vertices.col(2) + vertices.col(0));
      } else { // i == 2
        dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
      }
      BasisFunctionType weight =
          acc(lowestIndicesOfElementsAdjacentToEdges, edgeIndex) == elementIndex
              ? 1.
              : -1.;
      globalDofs.push_back(globalDofIndex);
      globalDofWeights.push_back(weight);
      if (globalDofIndex < 0)
        continue; // constrained edge
      acc(m_global2localDofs, globalDofIndex)
          .push_back(LocalDof(elementIndex, i));

      extendBoundingBox(acc(m_globalDofBoundingBoxes, globalDofIndex),
                        vertices);
      setBoundingBoxReference<CoordinateType>(
          acc(m_globalDofBoundingBoxes, globalDofIndex), dofPosition);

      ++flatLocalDofCount;
    }
    it->next();
  }

#ifndef NDEBUG
  for (size_t i = 0; i < m_globalDofBoundingBoxes.size(); ++i) {
    const BoundingBox<CoordinateType> &bbox = acc(m_globalDofBoundingBoxes, i);

    assert(bbox.reference.x >= bbox.lbound.x);
    assert(bbox.reference.y >= bbox.lbound.y);
    assert(bbox.reference.z >= bbox.lbound.z);
    assert(bbox.reference.x <= bbox.ubound.x);
    assert(bbox.reference.y <= bbox.ubound.y);
    assert(bbox.reference.z <= bbox.ubound.z);
  }
#endif // NDEBUG
*/
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      flatLocalDofCount, m_local2globalDofs, m_flatLocal2localDofs);

  //    // Iterate over elements
  //    std::unique_ptr<EntityIterator<elementCodim> > it =
  //        m_view->entityIterator<elementCodim>();
  //    m_flatLocalDofCount = 0;
  //    while (!it->finished())
  //    {
  //        const Entity<elementCodim>& element = it->entity();
  //        EntityIndex elementIndex = elementMapper.entityIndex(element);

  //        const int vertexCount = element.template
  // subEntityCount<vertexCodim>();
  //        const int edgeCount = vertexCount;
  //        m_flatLocalDofCount += edgeCount;

  //        // List of global DOF indices corresponding to the local DOFs of the
  //        // current element
  //        std::vector<GlobalDofIndex>& globalDofs =
  // m_local2globalDofs[elementIndex];
  //        // List of weights of the local DOFs residing on the current element
  //        std::vector<BasisFunctionType>& globalDofWeights =
  //            m_local2globalDofWeights[elementIndex];
  //        globalDofs.resize(edgeCount);
  //        globalDofWeights.resize(edgeCount);
  //        for (int i = 0; i < edgeCount; ++i)
  //        {
  //            GlobalDofIndex globalDofIndex =
  //                indexSet.subEntityIndex(element, i, edgeCodim);
  //            int vertex1Index, vertex2Index;
  //            // Handle Dune's funny subentity indexing order
  //            if (edgeCount == 3) {
  //                if (i == 0) {
  //                    vertex1Index = indexSet.subEntityIndex(element, 0,
  // vertexCodim);
  //                    vertex2Index = indexSet.subEntityIndex(element, 1,
  // vertexCodim);
  //                } else if (i == 1) {
  //                    vertex1Index = indexSet.subEntityIndex(element, 2,
  // vertexCodim);
  //                    vertex2Index = indexSet.subEntityIndex(element, 0,
  // vertexCodim);
  //                } else { // i == 2
  //                    vertex1Index = indexSet.subEntityIndex(element, 1,
  // vertexCodim);
  //                    vertex2Index = indexSet.subEntityIndex(element, 2,
  // vertexCodim);
  //                }
  //            } else // edgeCount == 4
  //                throw std::runtime_error(
  //                    "RaviartThomas0VectorSpaceBarycentric::"
  //                    "assignDofsImpl(): support for quadrilaterals not in
  // place yet");
  //            BasisFunctionType weight = vertex1Index < vertex2Index ? 1. :
  // -1.;
  //            globalDofs[i] = globalDofIndex;
  //            globalDofWeights[i] = weight;
  //            m_global2localDofs[globalDofIndex].push_back(LocalDof(elementIndex,
  // i));
  //            m_global2localDofWeights[globalDofIndex].push_back(weight);
  //        }
  //        it->next();
  //    }

  //    // Initialize the container mapping the flat local dof indices to
  //    // local dof indices
  //    m_flatLocal2localDofs.clear();
  //    m_flatLocal2localDofs.reserve(m_flatLocalDofCount);
  //    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
  //        for (size_t dof = 0; dof < m_local2globalDofs[e].size(); ++dof)
  //            m_flatLocal2localDofs.push_back(LocalDof(e, dof)); */
}
template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  const GridView &view = this->gridView();
  const Mapper &elementMapper = view.elementMapper();
  int index = elementMapper.entityIndex(element);
  if (m_elementIndex2Type[index] == Shapeset::TYPE1)
    return m_vectorBasisType1;
  else
    return m_vectorBasisType2;
}

template <typename BasisFunctionType>
size_t RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::globalDofCount() const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::flatLocalDofCount() const {
  return m_flatLocal2localDofs.size();
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getGlobalDofs(
    const Entity<0> &element, std::vector<GlobalDofIndex> &dofs,
    std::vector<BasisFunctionType> &dofWeights) const {
  const Mapper &mapper = m_view->elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = acc(m_local2globalDofs, index);
  dofWeights = acc(m_local2globalDofWeights, index);
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::global2localDofs(
    const std::vector<GlobalDofIndex> &globalDofs,
    std::vector<std::vector<LocalDof>> &localDofs,
    std::vector<std::vector<BasisFunctionType>> &localDofWeights) const {
  localDofs.resize(globalDofs.size());
  localDofWeights.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i) {
    acc(localDofs, i) = acc(m_global2localDofs, acc(globalDofs, i));
    std::vector<BasisFunctionType> &activeLdofWeights = acc(localDofWeights, i);
    activeLdofWeights.resize(localDofs[i].size());
    for (size_t j = 0; j < localDofs[i].size(); ++j) {
      LocalDof ldof = acc(localDofs[i], j);
      acc(activeLdofWeights, j) =
          acc(acc(m_local2globalDofWeights, ldof.entityIndex), ldof.dofIndex);
    }
  }
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::flatLocal2localDofs(
    const std::vector<FlatLocalDofIndex> &flatLocalDofs,
    std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    acc(localDofs, i) = acc(m_flatLocal2localDofs, acc(flatLocalDofs, i));
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getGlobalDofPositions(
    std::vector<Point3D<CoordinateType>> &positions) const {
  positions.resize(m_globalDofBoundingBoxes.size());
  for (size_t i = 0; i < m_globalDofBoundingBoxes.size(); ++i)
    acc(positions, i) = acc(m_globalDofBoundingBoxes, i).reference;
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getFlatLocalDofPositions(
    std::vector<Point3D<CoordinateType>> &positions) const {
  std::vector<BoundingBox<CoordinateType>> bboxes;
  getFlatLocalDofBoundingBoxes(bboxes);
  positions.resize(bboxes.size());
  for (size_t i = 0; i < bboxes.size(); ++i)
    acc(positions, i) = acc(bboxes, i).reference;
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getGlobalDofBoundingBoxes(
    std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  bboxes = m_globalDofBoundingBoxes;
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getFlatLocalDofBoundingBoxes(
    std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  const int flatLocalDofCount = m_flatLocal2localDofs.size();
  bboxes.resize(flatLocalDofCount);

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
          "RaviartThomas0VectorSpaceBarycentric::getFlatLocalDofBoundingBoxes(): "
          "only triangular elements are supported at present");
    it->next();
  }

  size_t flatLdofIndex = 0;
  Vector<CoordinateType> dofPosition;
  for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
    for (size_t v = 0; v < acc(m_local2globalDofs, e).size(); ++v)
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
  assert(flatLdofIndex == flatLocalDofCount);
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getGlobalDofNormals(
    std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = 2;
  const int worldDim = 3;
  const int globalDofCount_ = globalDofCount();
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

  for (size_t g = 0; g < globalDofCount_; ++g) {
    Point3D<CoordinateType> &normal = acc(normals, g);
    normal.x = 0.;
    normal.y = 0.;
    for (size_t l = 0; l < m_global2localDofs[g].size(); ++l) {
      normal.x += elementNormals(0, m_global2localDofs[g][l].entityIndex);
      normal.y += elementNormals(1, m_global2localDofs[g][l].entityIndex);
      normal.z += elementNormals(2, m_global2localDofs[g][l].entityIndex);
    }
    normal.x /= m_global2localDofs[g].size();
    normal.y /= m_global2localDofs[g].size();
    normal.z /= m_global2localDofs[g].size();
  }
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::getFlatLocalDofNormals(
    std::vector<Point3D<CoordinateType>> &normals) const {
  const int gridDim = 2;
  const int worldDim = 3;
  normals.resize(flatLocalDofCount());

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
      elementNormals(dim, index) = center(dim, 0);
    it->next();
  }

  size_t flatLdofIndex = 0;
  assert(m_local2globalDofs.size() == elementCount);
  for (size_t e = 0; e < elementCount; ++e)
    for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v)
      if (m_local2globalDofs[e][v] >= 0) { // is this LDOF used?
        normals[flatLdofIndex].x = elementNormals(0, e);
        normals[flatLdofIndex].y = elementNormals(1, e);
        normals[flatLdofIndex].z = elementNormals(2, e);
        ++flatLdofIndex;
      }
  assert(flatLdofIndex == flatLocalDofCount());
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::dumpClusterIds(
    const char *fileName,
    const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>::dumpClusterIdsEx(
    const char *fileName, const std::vector<unsigned int> &clusterIdsOfDofs,
    DofType dofType) const {
  throw std::runtime_error("RaviartThomas0VectorSpaceBarycentric::"
                           "dumpClusterIdsEx(): Not implemented yet");
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveRaviartThomas0VectorSpaceBarycentric(const shared_ptr<const Grid> &grid) {

  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType,
                        RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>>(grid));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveRaviartThomas0VectorSpaceBarycentric(const shared_ptr<const Grid> &grid,
                                     const std::vector<int> &domains,
                                     bool open) {

  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType,
                        RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>>(
          grid, domains, open));
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>>
adaptiveRaviartThomas0VectorSpaceBarycentric(const shared_ptr<const Grid> &grid,
                                     int domain, bool open) {

  return shared_ptr<Space<BasisFunctionType>>(
      new AdaptiveSpace<BasisFunctionType,
                        RaviartThomas0VectorSpaceBarycentric<BasisFunctionType>>(
          grid, std::vector<int>({domain}), open));
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)                                      \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveRaviartThomas0VectorSpaceBarycentric<BASIS>(const shared_ptr<const Grid> &); \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveRaviartThomas0VectorSpaceBarycentric<BASIS>(const shared_ptr<const Grid> &,  \
                                              const std::vector<int> &, bool); \
  template shared_ptr<Space<BASIS>>                                            \
  adaptiveRaviartThomas0VectorSpaceBarycentric<BASIS>(const shared_ptr<const Grid> &,  \
                                              int, bool)

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(RaviartThomas0VectorSpaceBarycentric);

} // namespace Bempp
