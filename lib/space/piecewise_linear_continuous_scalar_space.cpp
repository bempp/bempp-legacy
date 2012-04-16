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

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

#include <dune/localfunctions/lagrange/p1/p1localbasis.hh>
#include <dune/localfunctions/lagrange/q1/q1localbasis.hh>

#include <stdexcept>
#include <iostream>

#include "../grid/geometry.hpp"


namespace Bempp
{

template <typename ValueType>
PiecewiseLinearContinuousScalarSpace<ValueType>::
PiecewiseLinearContinuousScalarSpace(Grid& grid) :
     ScalarSpace<ValueType>(grid)
{
    const int gridDim = grid.dim();
    if (gridDim != 1 && gridDim != 2)
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "PiecewiseLinearContinuousScalarSpace(): "
                                    "only 1- and 2-dimensional grids are supported");
    m_view = this->m_grid.leafView();
}

template <typename ValueType>
int PiecewiseLinearContinuousScalarSpace<ValueType>::domainDimension() const
{
    return this->m_grid.dim();
}

template <typename ValueType>
int PiecewiseLinearContinuousScalarSpace<ValueType>::codomainDimension() const
{
    return 1;
}

template <typename ValueType>
void PiecewiseLinearContinuousScalarSpace<ValueType>::getBases(
        const std::vector<const EntityPointer<0>*>& elements,
        std::vector<const Fiber::Basis<ValueType>*>& bases) const
{
    const int elementCount = elements.size();
    bases.resize(elementCount);
    for (int i = 0; i < elementCount; ++i)
        switch (elementVariant(elements[i]->entity()))
        {
        case 3:
            bases[i] = &m_triangleBasis;
            break;
        case 4:
            bases[i] = &m_quadrilateralBasis;
            break;
        case 2:
            bases[i] = &m_lineBasis;
            break;
        }
}

template <typename ValueType>
const Fiber::Basis<ValueType>&
PiecewiseLinearContinuousScalarSpace<ValueType>::basis(
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
    }
}

template <typename ValueType>
ElementVariant PiecewiseLinearContinuousScalarSpace<ValueType>::elementVariant(
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

template <typename ValueType>
void PiecewiseLinearContinuousScalarSpace<ValueType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewiseLinearContinuousScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename ValueType>
void PiecewiseLinearContinuousScalarSpace<ValueType>::assignDofs()
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
//            arma::Col<ValueType> vertex;
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

template <typename ValueType>
bool PiecewiseLinearContinuousScalarSpace<ValueType>::dofsAssigned() const
{
    const int gridDim = domainDimension();
    return globalDofCount() == m_view->entityCount(gridDim);
}

template <typename ValueType>
int PiecewiseLinearContinuousScalarSpace<ValueType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename ValueType>
void PiecewiseLinearContinuousScalarSpace<ValueType>::globalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename ValueType>
void PiecewiseLinearContinuousScalarSpace<ValueType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (int i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename ValueType>
void PiecewiseLinearContinuousScalarSpace<ValueType>::globalDofPositions(
        std::vector<Point3D<ValueType> >& positions) const
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
            arma::Col<ValueType> vertex;
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
            arma::Col<ValueType> vertex;
            e.geometry().getCenter(vertex);

            positions[index].x = vertex(0);
            positions[index].y = vertex(1);
            positions[index].z = vertex(2);
            it->next();
        }
    }
}


#ifdef COMPILE_FOR_FLOAT
template class PiecewiseLinearContinuousScalarSpace<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class PiecewiseLinearContinuousScalarSpace<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class PiecewiseLinearContinuousScalarSpace<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class PiecewiseLinearContinuousScalarSpace<std::complex<double> >;
#endif

} // namespace Bempp
