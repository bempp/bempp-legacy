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
#include "../grid/vtk_writer.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
PiecewiseConstantScalarSpace<BasisFunctionType>::
PiecewiseConstantScalarSpace(Grid& grid) :
     ScalarSpace<BasisFunctionType>(grid)
{
}

template <typename BasisFunctionType>
int PiecewiseConstantScalarSpace<BasisFunctionType>::domainDimension() const
{
    return this->m_grid.dim();
}

template <typename BasisFunctionType>
int PiecewiseConstantScalarSpace<BasisFunctionType>::codomainDimension() const
{
    return 1;
}

template <typename BasisFunctionType>
const Fiber::Basis<BasisFunctionType>&
PiecewiseConstantScalarSpace<BasisFunctionType>::basis(
        const Entity<0>& element) const
{
    return m_basis;
}

template <typename BasisFunctionType>
ElementVariant PiecewiseConstantScalarSpace<BasisFunctionType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.isTriangle())
        return 3;
    else
        return 4;
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("PiecewiseConstantScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::assignDofs()
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

template <typename BasisFunctionType>
bool PiecewiseConstantScalarSpace<BasisFunctionType>::dofsAssigned() const
{
    return globalDofCount() == m_view->entityCount(0);
}

template <typename BasisFunctionType>
size_t PiecewiseConstantScalarSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t PiecewiseConstantScalarSpace<BasisFunctionType>::flatLocalDofCount() const
{
    return globalDofCount();
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::globalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    // Use the fact that each element contains exactly one DOF
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = LocalDof(flatLocalDofs[i], 0 /* local DOF #0 */);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::globalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = domainDimension();
    const int globalDofCount_ = globalDofCount();
    positions.resize(globalDofCount_);

    const Mapper& mapper = m_view->elementMapper();

    if (gridDim == 1)
        throw NotImplementedError("PiecewiseConstantScalarSpace::"
                                  "globalDofPositions(): "
                                  "not implemented for 2D yet.");
    else {
        std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
        while (!it->finished())
        {
            const Entity<0>& e = it->entity();
            int index = mapper.entityIndex(e);
            arma::Col<CoordinateType> center;
            e.geometry().getCenter(center);

            positions[index].x = center(0);
            positions[index].y = center(1);
            positions[index].z = center(2);
            it->next();
        }
    }
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::flatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    return globalDofPositions(positions);
}

template <typename BasisFunctionType>
void PiecewiseConstantScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIds) const
{
    const size_t idCount = clusterIds.size();
    if (idCount != globalDofCount())
        throw std::invalid_argument("PiecewiseConstantScalarSpace::"
                                    "dumpClusterIds(): incorrect dimension");

    std::auto_ptr<GridView> view = this->m_grid.leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();
    arma::Row<double> data(idCount);
    for (size_t i = 0; i < idCount; ++i)
        data(i) = clusterIds[i];
    vtkWriter->addCellData(data, "ids");
    vtkWriter->write(fileName);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);

} // namespace Bempp
