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

#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../grid/vtk_writer.hpp"

#include <dune/localfunctions/lagrange/p1/p1localbasis.hh>
#include <dune/localfunctions/lagrange/q1/q1localbasis.hh>

#include <stdexcept>
#include <iostream>

namespace Bempp
{

template <typename BasisFunctionType>
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
PiecewiseLinearContinuousScalarSpace(Grid& grid) :
     ScalarSpace<BasisFunctionType>(grid)
{
    const int gridDim = grid.dim();
    if (gridDim != 1 && gridDim != 2)
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "PiecewiseLinearContinuousScalarSpace(): "
                                    "only 1- and 2-dimensional grids are supported");
    m_view = this->m_grid.leafView();
}

template <typename BasisFunctionType>
int PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::domainDimension() const
{
    return this->m_grid.dim();
}

template <typename BasisFunctionType>
int PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::codomainDimension() const
{
    return 1;
}

template <typename BasisFunctionType>
const Fiber::Basis<BasisFunctionType>&
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::basis(
        const Entity<0>& element) const
{
    switch (elementVariant(element))
    {
    case 3:
        return m_triangleBasis;
    case 4:
        return m_quadrilateralBasis;
    case 2:
        return m_lineBasis;
    default:
        throw std::logic_error("PiecewiseLinearContinuousScalarSpace::basis(): "
                               "invalid element variant, this shouldn't happen!");
    }
}

template <typename BasisFunctionType>
ElementVariant PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::elementVariant(
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
        throw std::runtime_error("PiecewiseLinearContinuousScalarSpace::"
                                 "elementVariant(): invalid geometry type, "
                                 "this shouldn't happen!");
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewiseLinearContinuousScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::assignDofs()
{
    const int gridDim = domainDimension();

    const Mapper& elementMapper = m_view->elementMapper();
    const IndexSet& indexSet = m_view->indexSet();

    // Global DOF numbers will be identical with vertex indices.
    // Thus, the will be as many global DOFs as there are vertices.
    int globalDofCount_ = m_view->entityCount(this->m_grid.dim());
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
    int vertexCount;
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        EntityIndex elementIndex = elementMapper.entityIndex(element);

        if (gridDim == 1)
            vertexCount = element.template subEntityCount<1>();
        else // gridDim == 2
            vertexCount = element.template subEntityCount<2>();

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
}

template <typename BasisFunctionType>
bool PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::dofsAssigned() const
{
    const int gridDim = domainDimension();
    return globalDofCount() == m_view->entityCount(gridDim);
}

template <typename BasisFunctionType>
size_t PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::globalDofs(
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
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::globalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = domainDimension();
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
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIds) const
{
    const size_t idCount = clusterIds.size();
    if (idCount != globalDofCount())
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "dumpClusterIds(): incorrect dimension");

    std::auto_ptr<GridView> view = this->m_grid.leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();
    arma::Row<double> data(idCount);
    for (size_t i = 0; i < idCount; ++i)
        data(i) = clusterIds[i];
    vtkWriter->addVertexData(data, "ids");
    vtkWriter->write(fileName);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseLinearContinuousScalarSpace);

} // namespace Bempp
