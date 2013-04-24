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

#include "unit_scalar_space.hpp"

#include "../common/not_implemented_error.hpp"
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
UnitScalarSpace<BasisFunctionType>::
UnitScalarSpace(const shared_ptr<const Grid>& grid) :
     ScalarSpace<BasisFunctionType>(grid), m_view(grid->leafView())
{
    assignDofsImpl();
}

template <typename BasisFunctionType>
const Space<BasisFunctionType>&
UnitScalarSpace<BasisFunctionType>::discontinuousSpace() const
{
    throw std::runtime_error("UnitScalarSpace::discontinuousSpace(): "
                             "not implemented yet");
}

template <typename BasisFunctionType>
bool
UnitScalarSpace<BasisFunctionType>::isDiscontinuous() const
{
    return false;
}

template <typename BasisFunctionType>
int UnitScalarSpace<BasisFunctionType>::domainDimension() const
{
    return this->grid()->dim();
}

template <typename BasisFunctionType>
int UnitScalarSpace<BasisFunctionType>::codomainDimension() const
{
    return 1;
}

template <typename BasisFunctionType>
const Fiber::Basis<BasisFunctionType>&
UnitScalarSpace<BasisFunctionType>::basis(
        const Entity<0>& element) const
{
    return m_basis;
}

template <typename BasisFunctionType>
ElementVariant UnitScalarSpace<BasisFunctionType>::elementVariant(
        const Entity<0>& element) const
{
    GeometryType type = element.type();
    if (type.dim() == 1)
        return 2;
    if (type.isTriangle())
        return 3;
    else
        return 4;
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::setElementVariant(
        const Entity<0>& element, ElementVariant variant)
{
    if (variant != elementVariant(element))
        // for this space, the element variants are unmodifiable,
        throw std::runtime_error("UnitScalarSpace::"
                                 "setElementVariant(): invalid variant");
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::assignDofsImpl()
{
    const Mapper& mapper = m_view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();

    const int elementCount = m_view->entityCount(0);

    // List of local DOFs corresponding to a single global DOF.
    std::vector<LocalDof> localDofs;
    localDofs.reserve(elementCount);
    // List of global DOF indices corresponding to the local DOFs of a single
    // element
    std::vector<GlobalDofIndex> globalDofs(1);
    globalDofs[0] = 0;

    // (Re)initialise member variables
    m_local2globalDofs.clear();
    m_local2globalDofs.resize(elementCount);
    m_global2localDofs.clear();

    while (!it->finished()) {
        EntityIndex index = mapper.entityIndex(it->entity());
        localDofs.push_back(LocalDof(index, 0 /* local DOF #0 */));
        m_local2globalDofs[index] = globalDofs;
        it->next();
    }
    assert(localDofs.size() == elementCount);
    m_global2localDofs.push_back(localDofs);
}

template <typename BasisFunctionType>
size_t UnitScalarSpace<BasisFunctionType>::globalDofCount() const
{
    return m_global2localDofs.size();
}

template <typename BasisFunctionType>
size_t UnitScalarSpace<BasisFunctionType>::flatLocalDofCount() const
{
    return m_view->entityCount(0);
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::getGlobalDofs(
        const Entity<0>& element, std::vector<GlobalDofIndex>& dofs) const
{
    const Mapper& mapper = m_view->elementMapper();
    EntityIndex index = mapper.entityIndex(element);
    dofs = m_local2globalDofs[index];
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::global2localDofs(
        const std::vector<GlobalDofIndex>& globalDofs,
        std::vector<std::vector<LocalDof> >& localDofs) const
{
    localDofs.resize(globalDofs.size());
    for (size_t i = 0; i < globalDofs.size(); ++i)
        localDofs[i] = m_global2localDofs[globalDofs[i]];
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    // Use the fact that each element contains exactly one DOF
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = LocalDof(flatLocalDofs[i], 0 /* local DOF #0 */);
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    const int gridDim = domainDimension();
    assert(globalDofCount() == 1);
    positions.resize(1);
    positions[0].x = 0.;
    positions[0].y = 0.;
    positions[0].z = 0.;

    const Mapper& mapper = m_view->elementMapper();

    if (gridDim == 1)
        throw NotImplementedError(
                "UnitScalarSpace::globalDofPositions(): "
                "not implemented for 2D yet.");
    else {
        std::auto_ptr<EntityIterator<0> > it = m_view->entityIterator<0>();
        while (!it->finished())
        {
            const Entity<0>& e = it->entity();
            int index = mapper.entityIndex(e);
            arma::Col<CoordinateType> center;
            e.geometry().getCenter(center);

            positions[0].x += center(0);
            positions[0].y += center(1);
            positions[0].z += center(2);
            it->next();
        }
    }
    const int elementCount = m_view->entityCount(0);
    positions[0].x /= elementCount;
    positions[0].y /= elementCount;
    positions[0].z /= elementCount;
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    return getGlobalDofPositions(positions);
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfDofs) const
{
    dumpClusterIdsEx(fileName, clusterIdsOfDofs, GLOBAL_DOFS);
}

template <typename BasisFunctionType>
void UnitScalarSpace<BasisFunctionType>::dumpClusterIdsEx(
        const char* fileName,
        const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
        DofType dofType) const
{
    if (dofType != GLOBAL_DOFS && dofType != FLAT_LOCAL_DOFS)
        throw std::invalid_argument("PiecewiseConstantScalarSpace::"
                                    "dumpClusterIds(): invalid DOF type");
    const size_t idCount = clusterIdsOfGlobalDofs.size();
    if ((dofType == GLOBAL_DOFS && idCount != globalDofCount()) ||
            (dofType == FLAT_LOCAL_DOFS && idCount != flatLocalDofCount()))
        throw std::invalid_argument(
                "UnitScalarSpace::dumpClusterIds(): "
                "clusterIds has incorrect length");

    std::auto_ptr<GridView> view = this->grid()->leafView();
    std::auto_ptr<VtkWriter> vtkWriter = view->vtkWriter();
    arma::Row<double> data(idCount);
    for (size_t i = 0; i < idCount; ++i)
        data(i) = clusterIdsOfGlobalDofs[i];
    vtkWriter->addCellData(data, "ids");
    vtkWriter->write(fileName);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(UnitScalarSpace);

} // namespace Bempp
