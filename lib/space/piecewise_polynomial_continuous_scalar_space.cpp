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

#include "piecewise_polynomial_continuous_scalar_space.hpp"

#include "piecewise_polynomial_discontinuous_scalar_space.hpp"

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

#include <boost/array.hpp>

namespace Bempp
{

template <typename BasisFunctionType>
PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::
PiecewisePolynomialContinuousScalarSpace(const shared_ptr<const Grid>& grid,
                                         int polynomialOrder) :
    ScalarSpace<BasisFunctionType>(grid), m_polynomialOrder(polynomialOrder),
    m_flatLocalDofCount(0)
{
    const int gridDim = grid->dim();
    if (gridDim != 2)
        throw std::invalid_argument("PiecewisePolynomialContinuousScalarSpace::"
                                    "PiecewisePolynomialContinuousScalarSpace(): "
                                    "2-dimensional grids are supported");
    m_view = grid->leafView();
    if (polynomialOrder == 1)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 1>());
    else if (polynomialOrder == 2)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 2>());
    else if (polynomialOrder == 3)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 3>());
    else if (polynomialOrder == 4)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 4>());
    else if (polynomialOrder == 5)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 5>());
    else if (polynomialOrder == 6)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 6>());
    else if (polynomialOrder == 7)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 7>());
    else if (polynomialOrder == 8)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 8>());
    else if (polynomialOrder == 9)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 9>());
    else if (polynomialOrder == 10)
        m_triangleBasis.reset(new Fiber::LagrangeScalarBasis<3, BasisFunctionType, 10>());
    else
        throw std::invalid_argument("PiecewisePolynomialContinuousScalarSpace::"
                                    "PiecewisePolynomialContinuousScalarSpace(): "
                                    "polynomialOrder must be >= 1 and <= 10");
    assignDofsImpl();
}

template <typename BasisFunctionType>
PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::
~PiecewisePolynomialContinuousScalarSpace()
{
}

template <typename BasisFunctionType>
int PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::domainDimension() const
{
    return this->grid()->dim();
}

template <typename BasisFunctionType>
int PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::codomainDimension() const
{
    return 1;
}

template <typename BasisFunctionType>
const Fiber::Basis<BasisFunctionType>&
PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::basis(
        const Entity<0>& element) const
{
    if (elementVariant(element) == 3)
        return *m_triangleBasis;
    throw std::logic_error("PiecewisePolynomialContinuousScalarSpace::basis(): "
                           "invalid element variant, this shouldn't happen!");
}

template <typename BasisFunctionType>
ElementVariant PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.isLine())
        return 2;
    else if (type.isTriangle())
        return 3;
    else if (type.isQuadrilateral())
        return 4;
    else
        throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                                 "elementVariant(): invalid geometry type, "
                                 "this shouldn't happen!");
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
shared_ptr<const Space<BasisFunctionType> >
PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::discontinuousSpace(
    const shared_ptr<const Space<BasisFunctionType> >& self) const
{
   if (!m_discontinuousSpace) {
       tbb::mutex::scoped_lock lock(m_discontinuousSpaceMutex);
       typedef PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType>
               DiscontinuousSpace;
       if (!m_discontinuousSpace)
           m_discontinuousSpace.reset(new DiscontinuousSpace(
                                          this->grid(), m_polynomialOrder));
   }
   return m_discontinuousSpace;
}

template <typename BasisFunctionType>
bool
PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::isDiscontinuous() const
{
    return false;
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::assignDofsImpl()
{
    // TODO: refactor this function, it's way too long!

    // In addition to DOF assignment, this function also precalculates bounding 
    // boxes of global DOFs

    const int elementCount = m_view->entityCount(0);
    if (elementCount == 0)
        return;
    const int gridDim = this->domainDimension();
    if (gridDim != 2)
        throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                                 "assignDofsImpl(): only 2-dimensional grids "
                                 "are supported at present");
    const int vertexCodim = gridDim;
    const int edgeCodim = vertexCodim - 1;

    // const Mapper& elementMapper = m_view->elementMapper();
    const IndexSet& indexSet = m_view->indexSet();

    // Map vertices to global dofs
    const int vertexCount = m_view->entityCount(2);
    // For the time being all vertices are assumed unconstrained, so the array
    // below is unnecessary (it's just an identity mapping). But it will become
    // useful once we add support for mixed problems or triple junctions.
    std::vector<GlobalDofIndex> vertexGlobalDofs(vertexCount);
    for (int i = 0; i < vertexCount; ++i)
        acc(vertexGlobalDofs, (size_t)i) = i;
    m_vertexGlobalDofCount = vertexCount;

    // Map edges to global dofs
    const int edgeCount = m_view->entityCount(1);
    const int internalDofCountPerEdge = m_polynomialOrder - 1;
    // For the time being all edges are assumed unconstrained, so the array
    // below is unnecessary (it's just a linear mapping). But it will become
    // useful later.
    std::vector<GlobalDofIndex> edgeStartingGlobalDofs(edgeCount);
    for (int i = 0, gdof = m_vertexGlobalDofCount;
         i < edgeCount;
         gdof += internalDofCountPerEdge, ++i)
        acc(edgeStartingGlobalDofs, i) = gdof;
    m_edgeGlobalDofCount = edgeCount * internalDofCountPerEdge;

    // Map element interiors to global dofs
    const int bubbleDofCountPerTriangle =
        std::max(0, (m_polynomialOrder - 1) * (m_polynomialOrder - 2) / 2);
    const int bubbleDofCountPerQuad =
        std::max(0, (m_polynomialOrder - 1) * (m_polynomialOrder - 1));
    std::vector<int> bubbleDofCounts(elementCount);
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        EntityIndex elementIndex = indexSet.entityIndex(element);
        int vertexCount = element.template subEntityCount<2>();
        if (vertexCount != 3 && vertexCount != 4)
            throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                                     "assignDofsImpl(): elements must be "
                                     "triangular or quadrilateral");
        acc(bubbleDofCounts, elementIndex) =
            vertexCount == 3 ? bubbleDofCountPerTriangle : bubbleDofCountPerQuad;
        it->next();
    }
    std::vector<GlobalDofIndex> bubbleStartingGlobalDofs(elementCount);
    for (int i = 0, gdof = m_vertexGlobalDofCount + m_edgeGlobalDofCount;
         i < elementCount; gdof += acc(bubbleDofCounts, i), ++i)
        acc(bubbleStartingGlobalDofs, i) = gdof;
    m_bubbleGlobalDofCount =
        bubbleStartingGlobalDofs.back() + bubbleDofCounts.back() -
        (m_vertexGlobalDofCount + m_edgeGlobalDofCount);

    const int globalDofCount_ = m_vertexGlobalDofCount + m_edgeGlobalDofCount +
        m_bubbleGlobalDofCount;

    // Initialise DOF maps
    const int localDofCountPerTriangle =
        (m_polynomialOrder + 1) * (m_polynomialOrder + 2) / 2;
    const int localDofCountPerQuad =
        (m_polynomialOrder + 1) * (m_polynomialOrder + 1);
    m_local2globalDofs.clear();
    std::vector<GlobalDofIndex> prototypeGlobalDofs;
    prototypeGlobalDofs.reserve(localDofCountPerTriangle);
    m_local2globalDofs.resize(elementCount, prototypeGlobalDofs);
    m_global2localDofs.clear();
    // std::vector<LocalDof> prototypeLocalDofs;
    // prototypeLocalDofs.reserve(localDofCountPerTriangle);
    m_global2localDofs.resize(globalDofCount_/*, prototypeLocalDofs*/);

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
    it = m_view->entityIterator<0>();
    arma::Mat<CoordinateType> vertices;
    arma::Col<CoordinateType> dofPosition;
    m_flatLocalDofCount = 0;
    std::vector<int> gdofAccessCounts(globalDofCount_, 0);
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const Geometry& geo = element.geometry();
        EntityIndex elementIndex = indexSet.entityIndex(element);
        typedef arma::Col<CoordinateType> Col;

        geo.getCorners(vertices);
        int vertexCount = vertices.n_cols;

        // List of global DOF indices corresponding to the local DOFs of the
        // current element
        std::vector<GlobalDofIndex>& globalDofs = acc(m_local2globalDofs, elementIndex);
        if (vertexCount == 3) {
            std::vector<int> ldofAccessCounts(localDofCountPerTriangle, 0);
            boost::array<int, 3> vertexIndices;
            for (int i = 0; i < 3; ++i)
                acc(vertexIndices, i) = indexSet.subEntityIndex(element, i, vertexCodim);
            globalDofs.resize(localDofCountPerTriangle);
            // vertex dofs
            {
                int ldof, gdof;

                ldof = 0;
                gdof = vertexGlobalDofs[acc(vertexIndices, 0)];
                acc(globalDofs, ldof) = gdof;
                acc(m_global2localDofs, gdof).push_back(LocalDof(elementIndex, ldof));
                ++acc(ldofAccessCounts, ldof);
                ++acc(gdofAccessCounts, gdof);

                extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                setBoundingBoxReference<CoordinateType>(
                    acc(m_globalDofBoundingBoxes, gdof), vertices.col(0));

                ldof = m_polynomialOrder;
                gdof = vertexGlobalDofs[acc(vertexIndices, 1)];
                acc(globalDofs, ldof) = gdof;
                acc(m_global2localDofs, gdof).push_back(LocalDof(elementIndex, ldof));
                ++acc(ldofAccessCounts, ldof);
                ++acc(gdofAccessCounts, gdof);

                extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                setBoundingBoxReference<CoordinateType>(
                    acc(m_globalDofBoundingBoxes, gdof), vertices.col(1));

                ldof = localDofCountPerTriangle - 1;
                gdof = vertexGlobalDofs[acc(vertexIndices, 2)];
                acc(globalDofs, ldof) = gdof;
                acc(m_global2localDofs, gdof).push_back(LocalDof(elementIndex, ldof));
                ++acc(ldofAccessCounts, ldof);
                ++acc(gdofAccessCounts, gdof);

                extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                setBoundingBoxReference<CoordinateType>(
                    acc(m_globalDofBoundingBoxes, gdof), vertices.col(2));
            }

            // edge dofs
            if (m_polynomialOrder >= 2) {
                int start, end, step;
                int edgeIndex;

                edgeIndex = indexSet.subEntityIndex(element, 0, edgeCodim);
                dofPosition = 0.5 * (vertices.col(0) + vertices.col(1));
                if (acc(vertexIndices, 0) < acc(vertexIndices, 1)) {
                    start = acc(edgeStartingGlobalDofs, edgeIndex);
                    end = start + internalDofCountPerEdge;
                    step = 1;
                } else {
                    end = acc(edgeStartingGlobalDofs, edgeIndex) - 1;
                    start = end + internalDofCountPerEdge;
                    step = -1;
                }
                for (int ldof = 1, gdof = start; gdof != end; ++ldof, gdof += step) {
                    acc(globalDofs, ldof) = gdof;
                    acc(m_global2localDofs, gdof).push_back(LocalDof(elementIndex, ldof));
                    ++acc(ldofAccessCounts, ldof);
                    ++acc(gdofAccessCounts, gdof);

                    extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                    setBoundingBoxReference<CoordinateType>(
                        acc(m_globalDofBoundingBoxes, gdof), dofPosition);
                }

                edgeIndex = indexSet.subEntityIndex(element, 1, edgeCodim);
                dofPosition = 0.5 * (vertices.col(0) + vertices.col(2));
                if (acc(vertexIndices, 0) < acc(vertexIndices, 2)) {
                    start = acc(edgeStartingGlobalDofs, edgeIndex);
                    end = start + internalDofCountPerEdge;
                    step = 1;
                } else {
                    end = acc(edgeStartingGlobalDofs, edgeIndex) - 1;
                    start = end + internalDofCountPerEdge;
                    step = -1;
                }
                for (int ldofy = 1, gdof = start; gdof != end; ++ldofy, gdof += step) {
                    int ldof = ldofy * (m_polynomialOrder + 1) -
                        ldofy * (ldofy - 1) / 2;
                    acc(globalDofs, ldof) = gdof;
                    acc(m_global2localDofs, gdof).push_back(LocalDof(elementIndex, ldof));
                    ++acc(ldofAccessCounts, ldof);
                    ++acc(gdofAccessCounts, gdof);

                    extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                    setBoundingBoxReference<CoordinateType>(
                        acc(m_globalDofBoundingBoxes, gdof), dofPosition);
                }

                edgeIndex = indexSet.subEntityIndex(element, 2, edgeCodim);
                dofPosition = 0.5 * (vertices.col(1) + vertices.col(2));
                if (acc(vertexIndices, 1) < acc(vertexIndices, 2)) {
                    start = acc(edgeStartingGlobalDofs, edgeIndex);
                    end = start + internalDofCountPerEdge;
                    step = 1;
                } else {
                    end = acc(edgeStartingGlobalDofs, edgeIndex) - 1;
                    start = end + internalDofCountPerEdge;
                    step = -1;
                }
                for (int ldofy = 1, gdof = start; gdof != end; ++ldofy, gdof += step) {
                    int ldof = ldofy * (m_polynomialOrder + 1) -
                        ldofy * (ldofy - 1) / 2 + (m_polynomialOrder - ldofy);
                    acc(globalDofs, ldof) = gdof;
                    acc(m_global2localDofs, gdof).push_back(LocalDof(elementIndex, ldof));
                    ++acc(ldofAccessCounts, ldof);
                    ++acc(gdofAccessCounts, gdof);

                    extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                    setBoundingBoxReference<CoordinateType>(
                        acc(m_globalDofBoundingBoxes, gdof), dofPosition);
                }
            }

            // bubble dofs
            if (m_polynomialOrder >= 3) {
                dofPosition = (vertices.col(0) + vertices.col(1) +
                               vertices.col(2)) / 3.;
                for (int ldofy = 1, gdof = acc(bubbleStartingGlobalDofs, elementIndex);
                     ldofy < m_polynomialOrder; ++ldofy)
                    for (int ldofx = 1; ldofx + ldofy < m_polynomialOrder;
                         ++ldofx, ++gdof) {
                        int ldof = ldofy * (m_polynomialOrder + 1) -
                            ldofy * (ldofy - 1) / 2 + ldofx;
                        acc(globalDofs, ldof) = gdof;
                        acc(m_global2localDofs, gdof).push_back(
                            LocalDof(elementIndex, ldof));
                        ++acc(ldofAccessCounts, ldof);
                        ++acc(gdofAccessCounts, gdof);

                        extendBoundingBox(acc(m_globalDofBoundingBoxes, gdof), vertices);
                        setBoundingBoxReference<CoordinateType>(
                            acc(m_globalDofBoundingBoxes, gdof), dofPosition);
                    }
            }
            for (size_t i = 0; i < ldofAccessCounts.size(); ++i)
                assert(acc(ldofAccessCounts, i) == 1);
        } else
            throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                                     "assignDofsImpl(): quadrilateral elements "
                                     "are not supported yet");

        m_flatLocalDofCount += localDofCountPerTriangle;
        it->next();
    }
    // for (size_t i = 0; i < gdofAccessCounts.size(); ++i)
    //     std::cout << i << " " << acc(gdofAccessCounts, i) << "\n";

#ifndef NDEBUG

   for (size_t i = 0; i < globalDofCount_; ++i) {
       const BoundingBox<CoordinateType>& bbox = acc(m_globalDofBoundingBoxes, i);

       assert(bbox.reference.x >= bbox.lbound.x);
       assert(bbox.reference.y >= bbox.lbound.y);
       assert(bbox.reference.z >= bbox.lbound.z);
       assert(bbox.reference.x <= bbox.ubound.x);
       assert(bbox.reference.y <= bbox.ubound.y);
       assert(bbox.reference.z <= bbox.ubound.z);
   }
#endif // NDEBUG

    // Initialize the container mapping the flat local dof indices to
    // local dof indices
    m_flatLocal2localDofs.clear();
    m_flatLocal2localDofs.reserve(m_flatLocalDofCount);
    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
        for (size_t dof = 0; dof < acc(m_local2globalDofs, e).size(); ++dof)
            m_flatLocal2localDofs.push_back(LocalDof(e, dof));
}

template <typename BasisFunctionType>
size_t PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::flatLocalDofCount() const
{
    return m_flatLocalDofCount;
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::getGlobalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    positions.resize(m_globalDofBoundingBoxes.size());
    for (size_t i = 0; i < m_globalDofBoundingBoxes.size(); ++i)
        acc(positions, i) = acc(m_globalDofBoundingBoxes, i).reference;
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                             "getFlatLocalDofPositions(): not implemented yet");
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::getGlobalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
    bboxes = m_globalDofBoundingBoxes;
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::
getFlatLocalDofBoundingBoxes(
       std::vector<BoundingBox<CoordinateType> >& bboxes) const
{
    throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                             "getFlatLocalDofBoundingBoxes(): not implemented yet");
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::getGlobalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                             "getGlobalDofNormals(): not implemented yet");
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                             "getFlatLocalDofNormals(): not implemented yet");
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs) const
{
    dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void PiecewisePolynomialContinuousScalarSpace<BasisFunctionType>::dumpClusterIdsEx(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs,
        DofType dofType) const
{
    throw std::runtime_error("PiecewisePolynomialContinuousScalarSpace::"
                             "dumpClusterIdsEx(): not implemented yet");
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewisePolynomialContinuousScalarSpace);

} // namespace Bempp
