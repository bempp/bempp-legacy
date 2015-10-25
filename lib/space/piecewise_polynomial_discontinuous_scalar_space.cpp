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

#include "piecewise_polynomial_discontinuous_scalar_space.hpp"
#include "adaptive_space.hpp"

#include "space_helper.hpp"

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
#include "../grid/grid_segment.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <stdexcept>
#include <iostream>

#include <boost/array.hpp>

namespace Bempp {

namespace {

// Helper class that has a constructor without the space order
// (needed for adaptive space template)
template <typename BasisFunctionType, int order>
class PiecewisePolynomialDiscontinuousScalarSpaceHelper :
    public Bempp::PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType> {

public:

        PiecewisePolynomialDiscontinuousScalarSpaceHelper(const shared_ptr<const Grid>& grid, GridSegment segment):
                Bempp::PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>(grid,order, segment){}

        PiecewisePolynomialDiscontinuousScalarSpaceHelper(const shared_ptr<const Grid>& grid):
                Bempp::PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>(grid,order){}


    };

}

template <typename BasisFunctionType>
PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    PiecewisePolynomialDiscontinuousScalarSpace(
        const shared_ptr<const Grid> &grid, int polynomialOrder)
    : ScalarSpace<BasisFunctionType>(grid), m_polynomialOrder(polynomialOrder),
      m_flatLocalDofCount(0) {
  initialize(GridSegment::wholeGrid(*grid));
}

template <typename BasisFunctionType>
PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    PiecewisePolynomialDiscontinuousScalarSpace(
        const shared_ptr<const Grid> &grid, int polynomialOrder,
        const GridSegment &segment, int dofMode)
    : ScalarSpace<BasisFunctionType>(grid), m_polynomialOrder(polynomialOrder),
      m_flatLocalDofCount(0) {
  initialize(segment, dofMode);
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::initialize(
    const GridSegment &segment, int dofMode) {
  const int gridDim = this->grid()->dim();
  if (gridDim != 2)
    throw std::invalid_argument("PiecewisePolynomialDiscontinuousScalarSpace::"
                                "initialize(): "
                                "only 2-dimensional grids are supported");
  if (!(dofMode & (REFERENCE_POINT_ON_SEGMENT | ELEMENT_ON_SEGMENT)))
    throw std::invalid_argument("PiecewisePolynomialDiscontinuousScalarSpace::"
                                "initialize(): invalid dofMode");
  m_view = this->grid()->leafView();
  if (m_polynomialOrder == 0)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 0>());
  else if (m_polynomialOrder == 1)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 1>());
  else if (m_polynomialOrder == 2)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 2>());
  else if (m_polynomialOrder == 3)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 3>());
  else if (m_polynomialOrder == 4)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 4>());
  else if (m_polynomialOrder == 5)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 5>());
  else if (m_polynomialOrder == 6)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 6>());
  else if (m_polynomialOrder == 7)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 7>());
  else if (m_polynomialOrder == 8)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 8>());
  else if (m_polynomialOrder == 9)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 9>());
  else if (m_polynomialOrder == 10)
    m_triangleShapeset.reset(
        new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 10>());
  else
    throw std::invalid_argument(
        "PiecewisePolynomialDiscontinuousScalarSpace::"
        "PiecewisePolynomialDiscontinuousScalarSpace(): "
        "polynomialOrder must be >= 0 and <= 10");
  assignDofsImpl(segment, dofMode);
}

template <typename BasisFunctionType>
PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::~PiecewisePolynomialDiscontinuousScalarSpace() {}

template <typename BasisFunctionType>
int PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::domainDimension() const {
  return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::codomainDimension() const {
  return 1;
}

template <typename BasisFunctionType>
const Fiber::Shapeset<BasisFunctionType> &
PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::shapeset(
    const Entity<0> &element) const {
  if (elementVariant(element) == 3)
    return *m_triangleShapeset;
  throw std::logic_error(
      "PiecewisePolynomialDiscontinuousScalarSpace::shapeset(): "
      "invalid element variant, this shouldn't happen!");
}

template <typename BasisFunctionType>
bool PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::spaceIsCompatible(const Space<BasisFunctionType> &other)
    const {

  typedef PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>
      thisSpaceType;

  if (other.grid().get() != this->grid().get())
    return false;

  if (other.spaceIdentifier() == this->spaceIdentifier()) {
    // Try to typecast the other space down.
    const thisSpaceType &temp = dynamic_cast<const thisSpaceType &>(other);
    if (this->m_polynomialOrder == temp.m_polynomialOrder)
      return true;
    else
      return false;
  } else
    return false;
}

template <typename BasisFunctionType>
ElementVariant
PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::elementVariant(
    const Entity<0> &element) const {
  GeometryType type = element.type();
  if (type.isLine())
    return 2;
  else if (type.isTriangle())
    return 3;
  else if (type.isQuadrilateral())
    return 4;
  else
    throw std::runtime_error("PiecewisePolynomialDiscontinuousScalarSpace::"
                             "elementVariant(): invalid geometry type, "
                             "this shouldn't happen!");
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::setElementVariant(const Entity<0> &element,
                                          ElementVariant variant) {
  if (variant != elementVariant(element))
    // for this space, the element variants are unmodifiable,
    throw std::runtime_error("PiecewisePolynomialDiscontinuousScalarSpace::"
                             "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType>>
PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    discontinuousSpace(
        const shared_ptr<const Space<BasisFunctionType>> &self) const {
  if (self.get() != this)
    throw std::invalid_argument(
        "PiecewisePolynomialDiscontinuousScalarSpace::discontinuousSpace(): "
        "argument should be a shared pointer to *this");
  return self;
}

template <typename BasisFunctionType>
bool PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::isDiscontinuous() const {
  return true;
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::assignDofsImpl(const GridSegment &segment,
                                       int dofMode) {
  const int gridDim = this->domainDimension();
  if (gridDim != 2)
    throw std::runtime_error("PiecewisePolynomialDiscontinuousScalarSpace::"
                             "assignDofsImpl(): only 2-dimensional grids "
                             "are supported at present");
  const int vertexCodim = 2, edgeCodim = 1, elementCodim = 0;
  const Mapper &elementMapper = m_view->elementMapper();
  const IndexSet &indexSet = m_view->indexSet();

  int elementCount = m_view->entityCount(0);

  const int localDofCountPerTriangle =
      (m_polynomialOrder + 1) * (m_polynomialOrder + 2) / 2;
  const int localDofCountPerQuad =
      (m_polynomialOrder + 1) * (m_polynomialOrder + 1);

  BoundingBox<CoordinateType> model;
  model.lbound.x = std::numeric_limits<CoordinateType>::max();
  model.lbound.y = std::numeric_limits<CoordinateType>::max();
  model.lbound.z = std::numeric_limits<CoordinateType>::max();
  model.ubound.x = -std::numeric_limits<CoordinateType>::max();
  model.ubound.y = -std::numeric_limits<CoordinateType>::max();
  model.ubound.z = -std::numeric_limits<CoordinateType>::max();
  m_globalDofBoundingBoxes.reserve(localDofCountPerQuad * elementCount);

  // (Re)initialise DOF maps
  m_local2globalDofs.clear();
  m_local2globalDofs.resize(elementCount);
  m_global2localDofs.clear();
  // estimated number of global DOFs
  m_global2localDofs.reserve(localDofCountPerQuad * elementCount);

  // Fill in global<->local dof maps
  std::unique_ptr<EntityIterator<0>> it = m_view->entityIterator<0>();
  Matrix<CoordinateType> vertices;
  Vector<CoordinateType> dofPosition;
  GlobalDofIndex globalDofCount = 0;
  while (!it->finished()) {
    const Entity<0> &element = it->entity();
    EntityIndex elementIndex = elementMapper.entityIndex(element);
    bool elementContained = !(dofMode & ELEMENT_ON_SEGMENT) ||
                            segment.contains(elementCodim, elementIndex);
    typedef Vector<CoordinateType> Col;
    const Geometry &geo = element.geometry();
    geo.getCorners(vertices);
    int vertexCount = vertices.cols();
    int localDofCount =
        vertexCount == 3 ? localDofCountPerTriangle : localDofCountPerQuad;

    // List of global DOF indices corresponding to the local DOFs of the
    // current element
    std::vector<GlobalDofIndex> &globalDofs = m_local2globalDofs[elementIndex];
    globalDofs.resize(localDofCount);
    // GlobalDofIndex gdofStart = globalDofCount;
    // for (int i = 0; i < localDofCount; ++i) {
    //     globalDofs.push_back(globalDofCount);
    //     std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, i));
    //     m_global2localDofs.push_back(localDofs);
    //     ++globalDofCount;
    // }
    // GlobalDofIndex gdofEnd = globalDofCount;

    // Bounding boxes
    BoundingBox<CoordinateType> bbox = model;
    extendBoundingBox(bbox, vertices);
    // m_globalDofBoundingBoxes.insert(m_globalDofBoundingBoxes.end(),
    //                                 localDofCount, bbox);
    if (vertexCount == 3) {
      int subEntityIndex;
      int ldof;

      if (m_polynomialOrder == 0) {
        ldof = 0;
        if (segment.contains(0, elementIndex)) {
          acc(globalDofs, ldof) = globalDofCount;
          std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
          m_global2localDofs.push_back(localDofs);
          m_globalDofBoundingBoxes.push_back(bbox);
          setBoundingBoxReference<CoordinateType>(
              acc(m_globalDofBoundingBoxes, globalDofCount),
              (vertices.col(0) + vertices.col(1) + vertices.col(2)) / 3);
          ++globalDofCount;
        } else
          acc(globalDofs, ldof) = -1;
      } else {
        // vertex dofs
        ldof = 0;
        subEntityIndex = indexSet.subEntityIndex(element, 0, vertexCodim);
        if (elementContained &&
            (!(dofMode & REFERENCE_POINT_ON_SEGMENT) ||
             segment.contains(vertexCodim, subEntityIndex))) {
          acc(globalDofs, ldof) = globalDofCount;
          std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
          m_global2localDofs.push_back(localDofs);
          m_globalDofBoundingBoxes.push_back(bbox);
          setBoundingBoxReference<CoordinateType>(
              acc(m_globalDofBoundingBoxes, globalDofCount), vertices.col(0));
          ++globalDofCount;
        } else
          acc(globalDofs, ldof) = -1;

        ldof = m_polynomialOrder;
        subEntityIndex = indexSet.subEntityIndex(element, 1, vertexCodim);
        if (elementContained &&
            (!(dofMode & REFERENCE_POINT_ON_SEGMENT) ||
             segment.contains(vertexCodim, subEntityIndex))) {
          acc(globalDofs, ldof) = globalDofCount;
          std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
          m_global2localDofs.push_back(localDofs);
          m_globalDofBoundingBoxes.push_back(bbox);
          setBoundingBoxReference<CoordinateType>(
              acc(m_globalDofBoundingBoxes, globalDofCount), vertices.col(1));
          ++globalDofCount;
        } else
          acc(globalDofs, ldof) = -1;

        ldof = localDofCount - 1;
        subEntityIndex = indexSet.subEntityIndex(element, 2, vertexCodim);
        if (elementContained &&
            (!(dofMode & REFERENCE_POINT_ON_SEGMENT) ||
             segment.contains(vertexCodim, subEntityIndex))) {
          acc(globalDofs, ldof) = globalDofCount;
          std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
          m_global2localDofs.push_back(localDofs);
          m_globalDofBoundingBoxes.push_back(bbox);
          setBoundingBoxReference<CoordinateType>(
              acc(m_globalDofBoundingBoxes, globalDofCount), vertices.col(2));
          ++globalDofCount;
        } else
          acc(globalDofs, ldof) = -1;

        // edge dofs
        if (m_polynomialOrder >= 2) {

          subEntityIndex = indexSet.subEntityIndex(element, 0, edgeCodim);
          if (elementContained &&
              (!(dofMode & REFERENCE_POINT_ON_SEGMENT) ||
               segment.contains(edgeCodim, subEntityIndex))) {
            dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
            for (int ldof = 1; ldof < m_polynomialOrder; ++ldof) {
              acc(globalDofs, ldof) = globalDofCount;
              std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
              m_global2localDofs.push_back(localDofs);
              m_globalDofBoundingBoxes.push_back(bbox);
              setBoundingBoxReference<CoordinateType>(
                  acc(m_globalDofBoundingBoxes, globalDofCount), dofPosition);
              ++globalDofCount;
            }
          } else
            for (int ldof = 1; ldof < m_polynomialOrder; ++ldof) {
              acc(globalDofs, ldof) = -1;
            }

          subEntityIndex = indexSet.subEntityIndex(element, 1, edgeCodim);
          if (elementContained &&
              (!(dofMode & REFERENCE_POINT_ON_SEGMENT) ||
               segment.contains(edgeCodim, subEntityIndex))) {
            dofPosition = 0.5 * (vertices.col(0) + vertices.col(2));
            for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy) {
              int ldof =
                  ldofy * (m_polynomialOrder + 1) - ldofy * (ldofy - 1) / 2;
              acc(globalDofs, ldof) = globalDofCount;
              std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
              m_global2localDofs.push_back(localDofs);
              m_globalDofBoundingBoxes.push_back(bbox);
              setBoundingBoxReference<CoordinateType>(
                  acc(m_globalDofBoundingBoxes, globalDofCount), dofPosition);
              ++globalDofCount;
            }
          } else
            for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy) {
              int ldof =
                  ldofy * (m_polynomialOrder + 1) - ldofy * (ldofy - 1) / 2;
              acc(globalDofs, ldof) = -1;
            }

          subEntityIndex = indexSet.subEntityIndex(element, 2, edgeCodim);
          if (elementContained &&
              (!(dofMode & REFERENCE_POINT_ON_SEGMENT) ||
               segment.contains(edgeCodim, subEntityIndex))) {
            dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
            for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy) {
              int ldof = ldofy * (m_polynomialOrder + 1) -
                         ldofy * (ldofy - 1) / 2 + (m_polynomialOrder - ldofy);
              acc(globalDofs, ldof) = globalDofCount;
              std::vector<LocalDof> localDofs(1, LocalDof(elementIndex, ldof));
              m_global2localDofs.push_back(localDofs);
              m_globalDofBoundingBoxes.push_back(bbox);
              setBoundingBoxReference<CoordinateType>(
                  acc(m_globalDofBoundingBoxes, globalDofCount), dofPosition);
              ++globalDofCount;
            }
          } else
            for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy) {
              int ldof = ldofy * (m_polynomialOrder + 1) -
                         ldofy * (ldofy - 1) / 2 + (m_polynomialOrder - ldofy);
              acc(globalDofs, ldof) = -1;
            }

          // dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
          // for (int ldof = 1; ldof < m_polynomialOrder; ++ldof){
          //     setBoundingBoxReference<CoordinateType>(
          //         acc(m_globalDofBoundingBoxes, gdofStart + ldof),
          // dofPosition);
          // }
          // dofPosition = 0.5 * (vertices.col(0) + vertices.col(2));
          // for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy) {
          //     int ldof = ldofy * (m_polynomialOrder + 1) -
          //         ldofy * (ldofy - 1) / 2;
          //     setBoundingBoxReference<CoordinateType>(
          //         acc(m_globalDofBoundingBoxes, gdofStart + ldof),
          // dofPosition);
          // }
          // dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
          // for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy) {
          //     int ldof = ldofy * (m_polynomialOrder + 1) -
          //         ldofy * (ldofy - 1) / 2 + (m_polynomialOrder - ldofy);
          //     setBoundingBoxReference<CoordinateType>(
          //         acc(m_globalDofBoundingBoxes, gdofStart + ldof),
          // dofPosition);
          // }
        }
        // bubble dofs
        if (m_polynomialOrder >= 3) {
          if (segment.contains(elementCodim, elementIndex)) {
            dofPosition =
                (vertices.col(0) + vertices.col(1) + vertices.col(2)) / 3.;
            for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy)
              for (int ldofx = 1; ldofx + ldofy < m_polynomialOrder; ++ldofx) {
                int ldof = ldofy * (m_polynomialOrder + 1) -
                           ldofy * (ldofy - 1) / 2 + ldofx;
                acc(globalDofs, ldof) = globalDofCount;
                std::vector<LocalDof> localDofs(1,
                                                LocalDof(elementIndex, ldof));
                m_global2localDofs.push_back(localDofs);
                m_globalDofBoundingBoxes.push_back(bbox);
                setBoundingBoxReference<CoordinateType>(
                    acc(m_globalDofBoundingBoxes, globalDofCount), dofPosition);
                ++globalDofCount;
              }
          } else
            for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy)
              for (int ldofx = 1; ldofx + ldofy < m_polynomialOrder; ++ldofx) {
                int ldof = ldofy * (m_polynomialOrder + 1) -
                           ldofy * (ldofy - 1) / 2 + ldofx;
                acc(globalDofs, ldof) = -1;
              }

          // dofPosition = (vertices.col(0) + vertices.col(1) +
          //                vertices.col(2)) / 3.;
          // for (int ldofy = 1; ldofy < m_polynomialOrder; ++ldofy)
          //     for (int ldofx = 1; ldofx + ldofy < m_polynomialOrder; ++ldofx)
          // {
          //         int ldof = ldofy * (m_polynomialOrder + 1) -
          //             ldofy * (ldofy - 1) / 2 + ldofx;
          //         setBoundingBoxReference<CoordinateType>(
          //             acc(m_globalDofBoundingBoxes, gdofStart + ldof),
          //             dofPosition);
          //     }
        }
      }
    }

    it->next();
  }

  // Initialize the container mapping the flat local dof indices to
  // local dof indices
  SpaceHelper<BasisFunctionType>::initializeLocal2FlatLocalDofMap(
      m_flatLocalDofCount, m_local2globalDofs, m_flatLocal2localDofs);

  m_flatLocalDofCount = m_global2localDofs.size();

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
}

template <typename BasisFunctionType>
size_t
PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::globalDofCount()
    const {
  return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::flatLocalDofCount() const {
  return m_flatLocalDofCount;
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<
    BasisFunctionType>::getGlobalDofs(const Entity<0> &element,
                                      std::vector<GlobalDofIndex> &dofs) const {
  const Mapper &mapper = m_view->elementMapper();
  EntityIndex index = mapper.entityIndex(element);
  dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    global2localDofs(const std::vector<GlobalDofIndex> &globalDofs,
                     std::vector<std::vector<LocalDof>> &localDofs) const {
  localDofs.resize(globalDofs.size());
  for (size_t i = 0; i < globalDofs.size(); ++i)
    localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    flatLocal2localDofs(const std::vector<FlatLocalDofIndex> &flatLocalDofs,
                        std::vector<LocalDof> &localDofs) const {
  localDofs.resize(flatLocalDofs.size());
  for (size_t i = 0; i < flatLocalDofs.size(); ++i)
    localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  positions.resize(m_globalDofBoundingBoxes.size());
  for (size_t i = 0; i < m_globalDofBoundingBoxes.size(); ++i)
    acc(positions, i) = acc(m_globalDofBoundingBoxes, i).reference;
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType>> &positions) const {
  getGlobalDofPositions(positions);
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  bboxes = m_globalDofBoundingBoxes;
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofBoundingBoxes(
        std::vector<BoundingBox<CoordinateType>> &bboxes) const {
  getGlobalDofBoundingBoxes(bboxes);
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    getGlobalDofNormals(std::vector<Point3D<CoordinateType>> &normals) const {
  SpaceHelper<BasisFunctionType>::getGlobalDofNormals_defaultImplementation(
      *m_view, m_global2localDofs, normals);
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType>> &normals) const {
  getGlobalDofNormals(normals);
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    dumpClusterIds(const char *fileName,
                   const std::vector<unsigned int> &clusterIdsOfDofs) const {
  dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>::
    dumpClusterIdsEx(const char *fileName,
                     const std::vector<unsigned int> &clusterIdsOfDofs,
                     DofType dofType) const {
  throw std::runtime_error("PiecewisePolynomialDiscontinuousScalarSpace::"
                           "dumpClusterIdsEx(): not implemented yet");
}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewisePolynomialDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid,
        int order)
{

    if (order==0)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,0>>(grid));
    if (order==1)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,1>>(grid));
    if (order==2)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,2>>(grid));
    if (order==3)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,3>>(grid));
    if (order==4)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,4>>(grid));
    if (order==5)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,5>>(grid));
    if (order==6)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,6>>(grid));
    if (order==7)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,7>>(grid));
    if (order==8)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,8>>(grid));
    if (order==9)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,9>>(grid));
    if (order==10)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,10>>(grid));

    throw std::runtime_error("adaptivePiecewisePolynomialDiscontinuousScalarSpace(): Wrong order");

}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewisePolynomialDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid, int order,
        const std::vector<int>& domains, bool open)
{

    if (order==0)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,0>>(grid, domains, open));
    if (order==1)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,1>>(grid, domains, open));
    if (order==2)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,2>>(grid, domains, open));
    if (order==3)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,3>>(grid, domains, open));
    if (order==4)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,4>>(grid, domains, open));
    if (order==5)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,5>>(grid, domains, open));
    if (order==6)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,6>>(grid, domains, open));
    if (order==7)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,7>>(grid, domains, open));
    if (order==8)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,8>>(grid, domains, open));
    if (order==9)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,9>>(grid, domains, open));
    if (order==10)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,10>>(grid, domains, open));

    throw std::runtime_error("adaptivePiecewisePolynomialDiscontinuousScalarSpace(): Wrong order");

}

template <typename BasisFunctionType>
shared_ptr<Space<BasisFunctionType>> adaptivePiecewisePolynomialDiscontinuousScalarSpace(const shared_ptr<const Grid>& grid, int order, 
        int domain, bool open)
{
    if (order==0)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,0>>(grid, std::vector<int>({domain}), open));
    if (order==1)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,1>>(grid, std::vector<int>({domain}), open));
    if (order==2)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,2>>(grid, std::vector<int>({domain}), open));
    if (order==3)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,3>>(grid, std::vector<int>({domain}), open));
    if (order==4)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,4>>(grid, std::vector<int>({domain}), open));
    if (order==5)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,5>>(grid, std::vector<int>({domain}), open));
    if (order==6)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,6>>(grid, std::vector<int>({domain}), open));
    if (order==7)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,7>>(grid, std::vector<int>({domain}), open));
    if (order==8)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,8>>(grid, std::vector<int>({domain}), open));
    if (order==9)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,9>>(grid, std::vector<int>({domain}), open));
    if (order==10)
        return shared_ptr<Space<BasisFunctionType>>(
                 new AdaptiveSpace<BasisFunctionType, PiecewisePolynomialDiscontinuousScalarSpaceHelper<BasisFunctionType,10>>(grid, std::vector<int>({domain}), open));
    
    throw std::runtime_error("adaptivePiecewisePolynomialDiscontinuousScalarSpace(): Wrong order");
}

#define INSTANTIATE_FREE_FUNCTIONS(BASIS)   \
    template shared_ptr<Space<BASIS>> adaptivePiecewisePolynomialDiscontinuousScalarSpace<BASIS>( \
            const shared_ptr<const Grid>&, int); \
    template shared_ptr<Space<BASIS>> adaptivePiecewisePolynomialDiscontinuousScalarSpace<BASIS>( \
            const shared_ptr<const Grid>&, int, \
            const std::vector<int>&, bool); \
    template shared_ptr<Space<BASIS>> adaptivePiecewisePolynomialDiscontinuousScalarSpace<BASIS>( \
            const shared_ptr<const Grid>&, int, \
            int, bool) 


FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_FREE_FUNCTIONS);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    PiecewisePolynomialDiscontinuousScalarSpace);

} // namespace Bempp
