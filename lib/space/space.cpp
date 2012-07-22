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

//#include <boost/type_traits/is_same.hpp>
//#include <boost/static_assert.hpp>

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

    const Grid& grid = space.grid();
    std::auto_ptr<GridView> view = grid.leafView();
    const IndexSet& indexSet = view->indexSet();
    const size_t elementCount = view->entityCount(0);

    std::vector<std::vector<GlobalDofIndex> > gdofs;
    gdofs.resize(elementCount);
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& e = it->entity();
        int index = indexSet.entityIndex(e);
        space.globalDofs(e, gdofs[index]);
        it->next();
    }

    rows.clear();
    cols.clear();
    rows.reserve(ldofCount);
    cols.reserve(ldofCount);

    size_t flatLdofIndex = 0;
    for (size_t e = 0; e < gdofs.size(); ++e) {
        for (size_t v = 0; v < gdofs[e].size(); ++v) {
            rows.push_back(flatLdofIndex);
            cols.push_back(gdofs[e][v]);
            ++flatLdofIndex;
        }
    }
    assert(rows.size() == ldofCount);
    assert(cols.size() == ldofCount);

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
    Epetra_LocalMap rowMap(rowCount, 0 /* index_base */, comm);
    Epetra_LocalMap columnMap(columnCount, 0 /* index_base */, comm);
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
Space<BasisFunctionType>::Space(Grid& grid) :
    m_grid(grid)
{
}

template <typename BasisFunctionType>
Space<BasisFunctionType>::~Space()
{
}

template <typename BasisFunctionType>
void getAllBases(const Space<BasisFunctionType>& space,
        std::vector<const Fiber::Basis<BasisFunctionType>*>& bases)
{
    std::auto_ptr<GridView> view = space.grid().leafView();
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
shared_ptr<DiscreteBoundaryOperator<ResultType> >
constructOperatorMappingGlobalToFlatLocalDofs(const Space<BasisFunctionType>& space)
{
    shared_ptr<Epetra_CrsMatrix> mat =
            constructGlobalToFlatLocalDofsMappingEpetraMatrix(space);
    return boost::make_shared<DiscreteSparseBoundaryOperator<ResultType> >(
                mat);
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
constructOperatorMappingFlatLocalToGlobalDofs(const Space<BasisFunctionType>& space)
{
    shared_ptr<Epetra_CrsMatrix> mat =
            constructGlobalToFlatLocalDofsMappingEpetraMatrix(space);
    return boost::make_shared<DiscreteSparseBoundaryOperator<ResultType> >(
                mat, NO_SYMMETRY, TRANSPOSE);
}
#endif // WITH_TRILINOS

#define INSTANTIATE_getAllBases(BASIS) \
    template \
    void getAllBases(const Space<BASIS>& space, \
                std::vector<const Fiber::Basis<BASIS>*>& bases)

#define INSTANTIATE_constructOperators(BASIS, RESULT) \
    template \
    shared_ptr<DiscreteBoundaryOperator<RESULT> > \
            constructOperatorMappingGlobalToFlatLocalDofs(const Space<BASIS>& space); \
    template \
    shared_ptr<DiscreteBoundaryOperator<RESULT> > \
            constructOperatorMappingFlatLocalToGlobalDofs(const Space<BASIS>& space)

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Space);

FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_getAllBases);
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_constructOperators);

} // namespace Bempp
