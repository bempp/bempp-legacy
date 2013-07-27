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

#include "space.hpp"
#include "bempp/common/config_trilinos.hpp"

#include "../assembly/discrete_sparse_boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"

#ifdef WITH_TRILINOS
#include <Epetra_CrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>
#endif

namespace Bempp
{

namespace
{

#ifdef WITH_TRILINOS
template <typename BasisFunctionType>
void constructGlobalToFlatLocalDofsMappingVectors(
        const Space<BasisFunctionType>& space,
        std::vector<int>& rows, std::vector<int>& cols,
        std::vector<double>& values)
{
    const int ldofCount = space.flatLocalDofCount();

    std::auto_ptr<GridView> view = space.grid()->leafView();
    const IndexSet& indexSet = view->indexSet();
    const size_t elementCount = view->entityCount(0);

    std::vector<std::vector<GlobalDofIndex> > gdofs(elementCount);
    std::vector<std::vector<BasisFunctionType> > ldofWeights(elementCount);
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        space.getGlobalDofs(e, gdofs[index], ldofWeights[index]);
        it->next();
    }

    rows.clear();
    cols.clear();
    rows.reserve(ldofCount);
    cols.reserve(ldofCount);

    size_t flatLdofIndex = 0;
    for (size_t e = 0; e < gdofs.size(); ++e) {
        for (size_t v = 0; v < gdofs[e].size(); ++v) {
            int gdofIndex = gdofs[e][v];
            if (gdofIndex >= 0) {
                rows.push_back(flatLdofIndex);
                cols.push_back(gdofIndex);
                ++flatLdofIndex;
            }
        }
    }
    if (rows.size() != ldofCount || cols.size() != ldofCount)
        throw std::runtime_error(
                "constructGlobalToFlatLocalDofsMappingVectors(): "
                "internal error: the number of local DOFs is different from "
                "expected. Report this problem to BEM++ developers");

    std::vector<double> tmp(ldofCount, 1.);
    values.swap(tmp);
}

template <typename BasisFunctionType>
shared_ptr<Epetra_CrsMatrix>
constructGlobalToFlatLocalDofsMappingEpetraMatrix(
        const Space<BasisFunctionType>& space)
{
    std::vector<int> rows, cols;
    std::vector<double> values;
    constructGlobalToFlatLocalDofsMappingVectors(space, rows, cols, values);

    assert(rows.size() == cols.size());
    assert(cols.size() == values.size());
    const size_t entryCount = rows.size();

    const size_t rowCount = space.flatLocalDofCount();
    const size_t columnCount = space.globalDofCount();

    Epetra_SerialComm comm; // To be replaced once we begin to use MPI
    Epetra_LocalMap rowMap( (int) rowCount, 0 /* index_base */, comm);
    Epetra_LocalMap columnMap( (int) columnCount, 0 /* index_base */, comm);
    shared_ptr<Epetra_CrsMatrix> result = boost::make_shared<Epetra_CrsMatrix>(
                Copy, rowMap, columnMap, 1 /* entries per row */);

    for (size_t i = 0; i < entryCount; ++i) {
        int errorCode = result->InsertGlobalValues(
                    rows[i], 1 /* number of inserted entries */,
                    &values[i], &cols[i]);
        assert(errorCode == 0);
    }
    result->FillComplete(columnMap, rowMap);

    return result;
}
#endif // WITH_TRILINOS

} // namespace

template <typename BasisFunctionType>
Space<BasisFunctionType>::Space(const shared_ptr<const Grid>& grid) :
    m_grid(grid)
{
    if (!grid)
        throw std::invalid_argument("Space::Space(): grid must not be a null "
                                    "pointer");
}

template <typename BasisFunctionType>
Space<BasisFunctionType>::~Space()
{
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::assignDofs()
{
}

template <typename BasisFunctionType>
bool Space<BasisFunctionType>::dofsAssigned() const
{
    return true;
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::getGlobalDofs(
    const Entity<0>& element,
    std::vector<GlobalDofIndex>& dofs) const
{
    throw std::runtime_error("Space::getGlobalDofs(): the variant of this function "
                             "with two arguments is deprecated and not "
                             "implemented in new Space subclasses. Use the "
                             "three-argument variant instead");
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::getGlobalDofs(
    const Entity<0>& element,
    std::vector<GlobalDofIndex>& dofs,
    std::vector<BasisFunctionType>& localDofWeights) const
{
    getGlobalDofs(element, dofs);
    localDofWeights.resize(dofs.size(), 1.);
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs) const
{
    throw std::runtime_error("Space::global2localDofs(): the variant of this function "
                             "with two arguments is deprecated and not "
                             "implemented in new Space subclasses. Use the "
                             "three-argument variant instead");
}

template <typename BasisFunctionType>
void Space<BasisFunctionType>::global2localDofs(
            const std::vector<GlobalDofIndex>& globalDofs,
            std::vector<std::vector<LocalDof> >& localDofs,
            std::vector<std::vector<BasisFunctionType> >& localDofWeights) const
{
    global2localDofs(globalDofs, localDofs);
    localDofWeights.resize(localDofs.size());
    for (size_t igdof = 0; igdof < localDofs.size(); ++igdof)
        localDofWeights[igdof].resize(localDofs[igdof].size(),
                                      static_cast<BasisFunctionType>(1.));
}

template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType>& space,
        std::vector<const Fiber::Basis<BasisFunctionType>*>& bases)
{
    std::auto_ptr<GridView> view = space.grid()->leafView();
    const Mapper& mapper = view->elementMapper();
    const int elementCount = view->entityCount(0);

    bases.resize(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        bases[mapper.entityIndex(e)] = &space.basis(e);
        it->next();
    }
}

#ifdef WITH_TRILINOS
template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteSparseBoundaryOperator<ResultType> >
constructOperatorMappingGlobalToFlatLocalDofs(const Space<BasisFunctionType>& space)
{
    shared_ptr<Epetra_CrsMatrix> mat =
            constructGlobalToFlatLocalDofsMappingEpetraMatrix(space);
    return boost::make_shared<DiscreteSparseBoundaryOperator<ResultType> >(
                mat);
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteSparseBoundaryOperator<ResultType> >
constructOperatorMappingFlatLocalToGlobalDofs(const Space<BasisFunctionType>& space)
{
    shared_ptr<Epetra_CrsMatrix> mat =
            constructGlobalToFlatLocalDofsMappingEpetraMatrix(space);
    return boost::make_shared<DiscreteSparseBoundaryOperator<ResultType> >(
                mat, NO_SYMMETRY, TRANSPOSE);
}
#endif // WITH_TRILINOS

BEMPP_GCC_DIAG_OFF(deprecated-declarations);

template <typename BasisFunctionType>
void Space<BasisFunctionType>::dumpClusterIdsEx(
            const char* fileName,
            const std::vector<unsigned int>& clusterIdsOfGlobalDofs,
            DofType dofType) const
{
    if (dofType == GLOBAL_DOFS)
        return dumpClusterIds(fileName, clusterIdsOfGlobalDofs);
    else if (dofType == FLAT_LOCAL_DOFS)
        throw std::runtime_error("Space::dumpClusterIdsEx(): "
                                 "dumping of flat local DOF not supported");
    else
        throw std::invalid_argument("Space::dumpClusterIdsEx(): "
                                    "invalid DOF type");
}

BEMPP_GCC_DIAG_ON(deprecated-declarations);

#define INSTANTIATE_getAllBases(BASIS) \
    template \
    void getAllBases(const Space<BASIS>& space, \
                std::vector<const Fiber::Basis<BASIS>*>& bases)

#define INSTANTIATE_constructOperators(BASIS, RESULT) \
    template \
    shared_ptr<DiscreteSparseBoundaryOperator<RESULT> > \
            constructOperatorMappingGlobalToFlatLocalDofs(const Space<BASIS>& space); \
    template \
    shared_ptr<DiscreteSparseBoundaryOperator<RESULT> > \
            constructOperatorMappingFlatLocalToGlobalDofs(const Space<BASIS>& space)

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Space);

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_getAllBases);
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_constructOperators);

} // namespace Bempp
