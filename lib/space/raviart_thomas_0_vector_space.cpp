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

#include "raviart_thomas_0_vector_space.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
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

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct RaviartThomas0VectorSpace<BasisFunctionType>::Impl
{
    typedef Fiber::HdivFunctionValueFunctor<CoordinateType>
    TransformationFunctor;

    Impl() : transformations(TransformationFunctor())
    {}

    Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
};
/** \endcond */

template <typename BasisFunctionType>
RaviartThomas0VectorSpace<BasisFunctionType>::
RaviartThomas0VectorSpace(const shared_ptr<const Grid>& grid) :
    Base(grid), m_impl(new Impl), m_flatLocalDofCount(0)
{
    if (grid->dim() != 2)
        throw std::invalid_argument("RaviartThomas0VectorSpace::"
                                    "RaviartThomas0VectorSpace(): "
                                    "grid must be 2-dimensional");
    m_view = grid->leafView();
    assignDofsImpl();
}

template <typename BasisFunctionType>
RaviartThomas0VectorSpace<BasisFunctionType>::
~RaviartThomas0VectorSpace()
{
}

template <typename BasisFunctionType>
const typename RaviartThomas0VectorSpace<BasisFunctionType>::
CollectionOfBasisTransformations&
RaviartThomas0VectorSpace<BasisFunctionType>::shapeFunctionValue() const
{
    return m_impl->transformations;
}

template <typename BasisFunctionType>
int RaviartThomas0VectorSpace<BasisFunctionType>::domainDimension() const
{
    return 2;
}

template <typename BasisFunctionType>
int RaviartThomas0VectorSpace<BasisFunctionType>::codomainDimension() const
{
    return 3;
}

template <typename BasisFunctionType>
const Fiber::Basis<BasisFunctionType>&
RaviartThomas0VectorSpace<BasisFunctionType>::basis(
        const Entity<0>& element) const
{
    switch (elementVariant(element))
    {
    case 3:
        return m_triangleBasis;
    case 4:
        throw std::logic_error("RaviartThomas0VectorSpace::basis(): "
                               "quadrilateral elements are not supported yet");
    default:
        throw std::logic_error("RaviartThomas0VectorSpace::basis(): "
                               "invalid element variant, this shouldn't happen!");
    }
}

template <typename BasisFunctionType>
ElementVariant RaviartThomas0VectorSpace<BasisFunctionType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.isTriangle())
        return 3;
    else if (type.isQuadrilateral())
        return 4;
    else
        throw std::runtime_error("RaviartThomas0VectorSpace::"
                                 "elementVariant(): invalid geometry type, "
                                 "this shouldn't happen!");
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("RaviartThomas0VectorSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::assignDofsImpl()
{
    const int gridDim = domainDimension();

    const Mapper& elementMapper = m_view->elementMapper();
    const IndexSet& indexSet = m_view->indexSet();

    // Global DOF numbers will be identical with edge indices.
    // Thus, the will be as many global DOFs as there are edges.
    int globalDofCount_ = m_view->entityCount(1);
    int elementCount = m_view->entityCount(0);

    // (Re)initialise DOF maps
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_local2globalDofWeights.clear();
    m_local2globalDofWeights.resize(elementCount);
    m_global2localDofs.clear();
    m_global2localDofs.resize(globalDofCount_);
    m_global2localDofWeights.clear();
    m_global2localDofWeights.resize(globalDofCount_);
    // TODO: consider calling reserve(2) for each element of m_global2localDofs

    const int vertexCodim = 2;
    const int edgeCodim = 1;
    const int elementCodim = 0;

    // Iterate over elements
    std::auto_ptr<EntityIterator<elementCodim> > it =
        m_view->entityIterator<elementCodim>();
    m_flatLocalDofCount = 0;
    while (!it->finished())
    {
        const Entity<elementCodim>& element = it->entity();
        EntityIndex elementIndex = elementMapper.entityIndex(element);

        const int vertexCount = element.template subEntityCount<vertexCodim>();
        const int edgeCount = vertexCount;
        m_flatLocalDofCount += edgeCount;

        // List of global DOF indices corresponding to the local DOFs of the
        // current element
        std::vector<GlobalDofIndex>& globalDofs = m_local2globalDofs[elementIndex];
        // List of weights of the local DOFs residing on the current element
        std::vector<BasisFunctionType>& globalDofWeights =
            m_local2globalDofWeights[elementIndex];
        globalDofs.resize(edgeCount);
        globalDofWeights.resize(edgeCount);
        for (int i = 0; i < edgeCount; ++i)
        {
            GlobalDofIndex globalDofIndex =
                indexSet.subEntityIndex(element, i, edgeCodim);
            int vertex1Index, vertex2Index;
            // Handle Dune's funny subentity indexing order
            if (edgeCount == 3) {
                if (i == 0) {
                    vertex1Index = indexSet.subEntityIndex(element, 0, vertexCodim);
                    vertex2Index = indexSet.subEntityIndex(element, 1, vertexCodim);
                } else if (i == 1) {
                    vertex1Index = indexSet.subEntityIndex(element, 2, vertexCodim);
                    vertex2Index = indexSet.subEntityIndex(element, 0, vertexCodim);
                } else { // i == 2
                    vertex1Index = indexSet.subEntityIndex(element, 1, vertexCodim);
                    vertex2Index = indexSet.subEntityIndex(element, 2, vertexCodim);
                }
            } else // edgeCount == 4
                throw std::runtime_error(
                    "RaviartThomas0VectorSpace::"
                    "assignDofsImpl(): support for quadrilaterals not in place yet");
            BasisFunctionType weight = vertex1Index < vertex2Index ? 1. : -1.;
            globalDofs[i] = globalDofIndex;
            globalDofWeights[i] = weight;
            m_global2localDofs[globalDofIndex].push_back(LocalDof(elementIndex, i));
            m_global2localDofWeights[globalDofIndex].push_back(weight);
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
size_t RaviartThomas0VectorSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t RaviartThomas0VectorSpace<BasisFunctionType>::flatLocalDofCount() const
{
    return m_flatLocalDofCount;
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::getGlobalDofs(
        const Entity<0>& element,
        std::vector<GlobalDofIndex>& dofs,
        std::vector<BasisFunctionType>& dofWeights) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
    dofWeights = m_local2globalDofWeights[index];
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs,
        std::vector<std::vector<BasisFunctionType> >& localDofWeights) const
{
    localDofs.resize(globalDofs.size());
    localDofWeights.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofWeights[i] = m_global2localDofWeights[globalDofs[i]];
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = domainDimension();
    const int globalDofCount_ = globalDofCount();
    positions.resize(globalDofCount_);

    const IndexSet& indexSet = m_view->indexSet();

    std::auto_ptr<EntityIterator<1> > it = m_view->entityIterator<1>();
    while (!it->finished())
    {
        const Entity<1>& e = it->entity();
        int index = indexSet.entityIndex(e);
        arma::Col<CoordinateType> center;
        e.geometry().getCenter(center);

        positions[index].x = center(0);
        positions[index].y = center(1);
        positions[index].z = center(2);
        it->next();
    }
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = domainDimension();
    const int worldDim = this->grid()->dimWorld();
    positions.resize(m_flatLocalDofCount);

    const IndexSet& indexSet = m_view->indexSet();
    int elementCount = m_view->entityCount(0);

    arma::Mat<CoordinateType> elementCenters(worldDim, elementCount);
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
    arma::Col<CoordinateType> center;
    while (!it->finished())
    {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        e.geometry().getCenter(center);

        for (int dim = 0; dim < worldDim; ++dim)
            elementCenters(dim, index) = center(dim);
        it->next();
    }

    size_t flatLdofIndex = 0;
    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
        for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
            positions[flatLdofIndex].x = elementCenters(0, e);
            positions[flatLdofIndex].y = elementCenters(1, e);
            positions[flatLdofIndex].z = elementCenters(2, e);
            ++flatLdofIndex;
        }
    assert(flatLdofIndex == m_flatLocalDofCount);
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::getGlobalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    const int gridDim = domainDimension();
    const int worldDim = this->grid()->dimWorld();
    const int globalDofCount_ = globalDofCount();
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
void RaviartThomas0VectorSpace<BasisFunctionType>::getFlatLocalDofNormals(
        std::vector<Point3D<CoordinateType> >& normals) const
{
    const int gridDim = domainDimension();
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
            elementNormals(dim, index) = center(dim);
        it->next();
    }

    size_t flatLdofIndex = 0;
    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
        for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
            normals[flatLdofIndex].x = elementNormals(0, e);
            normals[flatLdofIndex].y = elementNormals(1, e);
            normals[flatLdofIndex].z = elementNormals(2, e);
            ++flatLdofIndex;
        }
    assert(flatLdofIndex == m_flatLocalDofCount);
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs) const
{
    dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void RaviartThomas0VectorSpace<BasisFunctionType>::dumpClusterIdsEx(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs,
        DofType dofType) const
{
    throw std::runtime_error("RaviartThomas0VectorSpace::"
                             "dumpClusterIdsEx(): Not implemented yet");
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(RaviartThomas0VectorSpace);

} // namespace Bempp
