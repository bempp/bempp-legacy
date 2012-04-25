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

#include "piecewise_constant_scalar_space.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"


namespace Bempp
{

template <typename ValueType>
PiecewiseConstantScalarSpace<ValueType>::PiecewiseConstantScalarSpace(Grid& grid) :
     ScalarSpace<ValueType>(grid)
{
}

template <typename ValueType>
int PiecewiseConstantScalarSpace<ValueType>::domainDimension() const
{
    return this->m_grid.dim();
}

template <typename ValueType>
int PiecewiseConstantScalarSpace<ValueType>::codomainDimension() const
{
    return 1;
}

template <typename ValueType>
void PiecewiseConstantScalarSpace<ValueType>::getBases(
        const std::vector<const EntityPointer<0>*>& elements,
        std::vector<const Fiber::Basis<ValueType>*>& bases) const
{
    const int elementCount = elements.size();
    bases.resize(elementCount);
    for (int i = 0; i < elementCount; ++i)
        bases[i] = &m_basis;
}

template <typename ValueType>
const Fiber::Basis<ValueType>& PiecewiseConstantScalarSpace<ValueType>::basis(
        const Entity<0>& element) const
{
    return m_basis;
}

//template <typename ValueType>
//void PiecewiseConstantScalarSpace<ValueType>::evaluateBasisFunctions(
//        ElementVariant elementVariant,
//        const arma::Mat<ctype>& local,
//        arma::Cube<ValueType>& result) const
//{
//    const int gridDim = domainDimension();
//    const int pointCount = local.n_cols;
//    const int codomainDim = codomainDimension();
//    const int functionCount = basisFunctionCount(elementVariant);

//    if (local.n_rows != gridDim)
//        throw std::invalid_argument(
//                "PiecewiseLinearContinuousScalarSpace::evaluateBasisFunctions(): "
//                "coordinates in 'local' have an incorrect number of components");

//    result.set_size(codomainDim, functionCount, pointCount);
//    // Basis function equal to 1. everywhere
//    result.fill(1.);
//}

//template <typename ValueType>
//void PiecewiseConstantScalarSpace<ValueType>::evaluateBasisFunctionDerivative(
//        ElementVariant elementVariant,
//        const arma::Mat<ctype>& local,
//        int direction,
//        arma::Cube<ValueType>& result) const
//{
//    const int gridDim = domainDimension();
//    const int pointCount = local.n_cols;
//    const int codomainDim = codomainDimension();
//    const int functionCount = basisFunctionCount(elementVariant);

//    assert(local.n_rows == gridDim);
//    assert(0 <= direction && direction < gridDim);

//    result.set_size(codomainDim, functionCount, pointCount);
//    // Constant function -> gradient null everywhere
//    result.fill(0.);
//}

template <typename ValueType>
ElementVariant PiecewiseConstantScalarSpace<ValueType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.isTriangle())
        return 3;
    else
        return 4;
}

template <typename ValueType>
void PiecewiseConstantScalarSpace<ValueType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewiseConstantScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename ValueType>
void PiecewiseConstantScalarSpace<ValueType>::assignDofs()
{
    m_view = this->m_grid.leafView();
    const Mapper& mapper = m_view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();

    // List of local DOFs corresponding to a single global DOF.
    std::vector<LocalDof> localDofs(1);
    // List of global DOF indices corresponding to the local DOFs of a single
    // element
    std::vector<GlobalDofIndex> globalDofs(1);
    // For this space, there is a one-to-one mapping between the local and
    // global DOFs, so the above vectors consist of one element only.

    // (Re)initialise member variables
    const int elementCount = m_view->entityCount(0);
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_global2localDofs.clear();
    m_global2localDofs.reserve(elementCount);

    int globalDofCount_ = 0;
    while (!it->finished())
    {
        EntityIndex index = mapper.entityIndex(it->entity());
        globalDofs[0] = globalDofCount_++;
        localDofs[0] = LocalDof(index, 0 /* local DOF #0 */);
        m_local2globalDofs[index] = globalDofs;
        m_global2localDofs.push_back(localDofs);
        it->next();
    }
    assert(globalDofCount_ == elementCount);
}

template <typename ValueType>
bool PiecewiseConstantScalarSpace<ValueType>::dofsAssigned() const
{
    return globalDofCount() == m_view->entityCount(0);
}

template <typename ValueType>
int PiecewiseConstantScalarSpace<ValueType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename ValueType>
void PiecewiseConstantScalarSpace<ValueType>::globalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename ValueType>
void PiecewiseConstantScalarSpace<ValueType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (int i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename ValueType>
void PiecewiseConstantScalarSpace<ValueType>::globalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = domainDimension();
    const int globalDofCount_ = globalDofCount();
    positions.resize(globalDofCount_);

    const IndexSet& indexSet = m_view->indexSet();

    if (gridDim == 1)
        throw NotImplementedError("PiecewiseConstantScalarSpace::"
                                  "globalDofPositions(): "
                                  "not implemented for 2D yet.");
    else {
        std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
        while (!it->finished())
        {
            const Entity<0>& e = it->entity();
            int index = indexSet.entityIndex(e);
            arma::Col<CoordinateType> center;
            e.geometry().getCenter(center);

            positions[index].x = center(0);
            positions[index].y = center(1);
            positions[index].z = center(2);
            it->next();
        }
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);

} // namespace Bempp
