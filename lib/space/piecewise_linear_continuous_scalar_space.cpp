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

#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
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
    ScalarSpace<BasisFunctionType>(grid), m_flatLocalDofCount(0)
{
    const int gridDim = grid.dim();
    if (gridDim != 1 && gridDim != 2)
        throw std::invalid_argument("PiecewiseLinearContinuousScalarSpace::"
                                    "PiecewiseLinearContinuousScalarSpace(): "
                                    "only 1- and 2-dimensional grids are supported");
    m_view = this->m_grid.leafView();
}

template <typename BasisFunctionType>
PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
~PiecewiseLinearContinuousScalarSpace()
{
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
    m_flatLocalDofCount = 0;
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        EntityIndex elementIndex = elementMapper.entityIndex(element);

        int vertexCount;
        if (gridDim == 1)
            vertexCount = element.template subEntityCount<1>();
        else // gridDim == 2
            vertexCount = element.template subEntityCount<2>();
        m_flatLocalDofCount += vertexCount;

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

    // Initialize the container mapping the flat local dof indices to 
    // local dof indices
    m_flatLocal2localDofs.clear();
    m_flatLocal2localDofs.reserve(m_flatLocalDofCount);
    for (size_t e = 0; e < m_local2globalDofs.size(); ++e)
        for (size_t dof = 0; dof < m_local2globalDofs[e].size(); ++dof)
            m_flatLocal2localDofs.push_back(LocalDof(e, dof));
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
size_t PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::flatLocalDofCount() const
{
// This is the correct implementation. Include it once the
// bug in FoamGrid is fixes.
//    if (gridDim == 1)
//        return m_view->entityCount(0) * 2;
//    else // gridDim == 2
//        return m_view->entityCount(GeometryType(GeometryType::cube, 2)) * 4 +
//                m_view->entityCount(GeometryType(GeometryType::simplex, 2)) * 3;

    return m_flatLocalDofCount;
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getGlobalDofs(
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
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::flatLocal2localDofs(
        const std::vector<FlatLocalDofIndex>& flatLocalDofs,
        std::vector<LocalDof>& localDofs) const
{
    localDofs.resize(flatLocalDofs.size());
    for (size_t i = 0; i < flatLocalDofs.size(); ++i)
        localDofs[i] = m_flatLocal2localDofs[flatLocalDofs[i]];
}

//template <typename BasisFunctionType>
//shared_ptr<DiscreteBoundaryOperator<
//typename PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::CoordinateType> >
//PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
//global2localDofsOperatorRealImpl() const
//{
//    std::vector<int> rows, cols;
//    std::vector<double> values;
//    constructGlobal2localDofsMappingVectors(rows, cols, values);
//    return boost::make_shared<DiscreteSparseBoundaryOperator<CoordinateType> >(
//                rows, cols, values);
//}

//template <typename BasisFunctionType>
//shared_ptr<DiscreteBoundaryOperator<
//typename PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::ComplexType> >
//PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
//global2localDofsOperatorComplexImpl() const
//{
//    std::vector<int> rows, cols;
//    std::vector<double> values;
//    constructGlobal2localDofsMappingVectors(rows, cols, values);
//    return boost::make_shared<DiscreteSparseBoundaryOperator<ComplexType> >(
//                rows, cols, values);
//}

//template <typename BasisFunctionType>
//shared_ptr<DiscreteBoundaryOperator<
//typename PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::CoordinateType> >
//PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
//local2globalDofsOperatorRealImpl() const
//{
//    std::vector<int> rows, cols;
//    std::vector<double> values;
//    constructGlobal2localDofsMappingVectors(rows, cols, values);
//    return boost::make_shared<DiscreteSparseBoundaryOperator<CoordinateType> >(
//                rows, cols, values, NO_SYMMETRY, TRANSPOSE);
//}

//template <typename BasisFunctionType>
//shared_ptr<DiscreteBoundaryOperator<
//typename PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::ComplexType> >
//PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
//local2globalDofsOperatorRealImpl() const
//{
//    std::vector<int> rows, cols;
//    std::vector<double> values;
//    constructGlobal2localDofsMappingVectors(rows, cols, values);
//    return boost::make_shared<DiscreteSparseBoundaryOperator<ComplexType> >(
//                rows, cols, values, NO_SYMMETRY, TRANSPOSE);
//}

//template <typename BasisFunctionType>
//void
//PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::
//constructGlobal2localDofsMappingVectors(
//        std::vector<int>& rows, std::vector<int>& cols,
//        std::vector<double>& values) const
//{
//    const size_t ldofCount = m_flatLocalDofCount;
//    rows.clear();
//    cols.clear();
//    rows.reserve(ldofCount);
//    cols.reserve(ldofCount);

//    size_t flatLdofIndex = 0;
//    for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
//        for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
//            rows.push_back(flatLdofIndex);
//            cols.push_back(m_local2globalDofs[e][v]);
//            ++flatLdofIndex;
//        }
//    }
//    assert(rows.size() == ldofCount);
//    assert(cols.size() == ldofCount);

//    std::vector<double> tmp(ldofCount, 1.);
//    values.swap(tmp);
//}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getGlobalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    if (!dofsAssigned())
        throw std::runtime_error(
                "PiecewiseLinearContinuousScalarSpace::getGlobalDofPositions(): "
                "assignDofs() must be called before calling this function");

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
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::getFlatLocalDofPositions(
        std::vector<Point3D<CoordinateType> >& positions) const
{
    if (!dofsAssigned())
        throw std::runtime_error(
                "PiecewiseLinearContinuousScalarSpace::getFlatLocalDofPositions(): "
                "assignDofs() must be called before calling this function");

    const int gridDim = domainDimension();
    const int worldDim = this->m_grid.dimWorld();
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
    if (gridDim == 1)
        for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
            for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
                positions[flatLdofIndex].x = elementCenters(0, e);
                positions[flatLdofIndex].y = elementCenters(1, e);
                positions[flatLdofIndex].z = 0.;
                ++flatLdofIndex;
            }
        }
    else // gridDim == 2
        for (size_t e = 0; e < m_local2globalDofs.size(); ++e) {
            for (size_t v = 0; v < m_local2globalDofs[e].size(); ++v) {
                positions[flatLdofIndex].x = elementCenters(0, e);
                positions[flatLdofIndex].y = elementCenters(1, e);
                positions[flatLdofIndex].z = elementCenters(2, e);
                ++flatLdofIndex;
            }
        }
    assert(flatLdofIndex == m_flatLocalDofCount);
}

template <typename BasisFunctionType>
void PiecewiseLinearContinuousScalarSpace<BasisFunctionType>::dumpClusterIds(
        const char* fileName,
        const std::vector<unsigned int>& clusterIds) const
{
    if (!dofsAssigned())
        throw std::runtime_error(
                "PiecewiseLinearContinuousScalarSpace::dumpClusterIds(): "
                "assignDofs() must be called before calling this function");

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
