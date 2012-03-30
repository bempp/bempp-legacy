#include "identity_operator.hpp"

#include <stdexcept>
#include <vector>

#include "assembly_options.hpp"
#include "discrete_dense_scalar_valued_linear_operator.hpp"
#include "discrete_sparse_scalar_valued_linear_operator.hpp"

#include "../common/types.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#ifdef WITH_TRILINOS
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>
#endif // WITH_TRILINOS

namespace Bempp
{

#ifdef WITH_TRILINOS
// Internal helper functions for Epetra
namespace
{

template <typename ValueType>
inline int epetraSumIntoGlobalValues(Epetra_FECrsMatrix& matrix,
                              const std::vector<int>& rowIndices,
                              const std::vector<int>& colIndices,
                              const arma::Mat<ValueType>& values)
{
    assert(rowIndices.size() == values.n_rows);
    assert(colIndices.size() == values.n_cols);
    // Convert data in ValueType into double (expected by Epetra)
    arma::Mat<double> doubleValues(values.n_rows, values.n_cols);
    for (int col = 0; col < values.n_cols; ++col)
        for (int row = 0; row < values.n_rows; ++row)
            doubleValues(row, col) = values(row, col);
    return matrix.SumIntoGlobalValues(rowIndices.size(), &rowIndices[0],
                                      colIndices.size(), &colIndices[0],
                                      doubleValues.memptr(),
                                      Epetra_FECrsMatrix::COLUMN_MAJOR);
}

// Specialisation for double -- no intermediate array is needed
template <>
inline int epetraSumIntoGlobalValues<double>(Epetra_FECrsMatrix& matrix,
                                      const std::vector<int>& rowIndices,
                                      const std::vector<int>& colIndices,
                                      const arma::Mat<double>& values)
{
    assert(rowIndices.size() == values.n_rows);
    assert(colIndices.size() == values.n_cols);
    return matrix.SumIntoGlobalValues(rowIndices.size(), &rowIndices[0],
                                      colIndices.size(), &colIndices[0],
                                      values.memptr(),
                                      Epetra_FECrsMatrix::COLUMN_MAJOR);
}

// TODO: decide whether/how to handle ValueType = complex<...>

} // anonymous namespace
#endif

template <typename ValueType>
bool IdentityOperator<ValueType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (repr == AssemblyOptions::DENSE || repr == AssemblyOptions::SPARSE);
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
IdentityOperator<ValueType>::assembleOperator(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        const typename IdentityOperator<ValueType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    throw std::runtime_error("IdentityOperator::assembleOperator(): "
                             "not implemented yet");
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
IdentityOperator<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const typename IdentityOperator<ValueType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    if (!testSpace.dofsAssigned() || !trialSpace.dofsAssigned())
        throw std::runtime_error("IdentityOperator::assembleWeakForm(): "
                                 "degrees of freedom must be assigned "
                                 "before calling assembleWeakForm()");
    if (&testSpace.grid() != &trialSpace.grid())
        throw std::runtime_error("IdentityOperator::assembleWeakForm(): "
                                 "testSpace and trialSpace must be defined over "
                                 "the same grid");

    // Prepare local assembler

    std::auto_ptr<GridView> view = trialSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Gather geometric data
    Fiber::RawGridGeometry<ValueType> rawGeometry;
    view->getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData());

    // Make geometry factory
    std::auto_ptr<GeometryFactory> geometryFactory =
            testSpace.grid().elementGeometryFactory();

    // Get pointers to test and trial bases of each element
    std::vector<const Fiber::Basis<ValueType>*> testBases;
    std::vector<const Fiber::Basis<ValueType>*> trialBases;
    testBases.reserve(elementCount);
    trialBases.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        testBases.push_back(&testSpace.basis(element));
        trialBases.push_back(&trialSpace.basis(element));
        it->next();
    }

    Fiber::OpenClHandler<ValueType,int> openClHandler(options.openClOptions());
    if (openClHandler.UseOpenCl())
        openClHandler.pushGeometry (rawGeometry.vertices(),
				    rawGeometry.elementCornerIndices());

    // Now create the assembler
    std::auto_ptr<LocalAssembler> assembler =
            factory.make(*geometryFactory, rawGeometry,
                         testBases, trialBases,
                         m_expression, m_expression, this->multiplier(),
                         openClHandler);

    return assembleWeakFormInternal(testSpace, trialSpace, *assembler, options);
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
IdentityOperator<ValueType>::assembleWeakFormInternal(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    switch (options.operatorRepresentation())
    {
    case AssemblyOptions::DENSE:
        return assembleWeakFormInDenseMode(
                    testSpace, trialSpace, assembler, options);
    case AssemblyOptions::SPARSE:
        return assembleWeakFormInSparseMode(
                    testSpace, trialSpace, assembler, options);
    default:
        throw std::runtime_error("IdentityOperator::assembleWeakForm(): "
                                 "invalid assembly mode");
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
IdentityOperator<ValueType>::assembleWeakFormInDenseMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        typename IdentityOperator<ValueType>::LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    // Fill local submatrices
    std::auto_ptr<GridView> view = testSpace.grid().leafView();
    const int elementCount = view->entityCount(0);
    std::vector<int> elementIndices(elementCount);
    for (int i = 0; i < elementCount; ++i)
        elementIndices[i] = i;
    std::vector<arma::Mat<ValueType> > localResult;
    assembler.evaluateLocalWeakForms(elementIndices, localResult);

    // Create the operator's matrix
    arma::Mat<ValueType> result(testSpace.globalDofCount(),
                                trialSpace.globalDofCount());
    result.fill(0.);

    // Retrieve global DOFs corresponding to local DOFs on all elements
    std::vector<std::vector<GlobalDofIndex> > trialGdofs(elementCount);
    std::vector<std::vector<GlobalDofIndex> > testGdofs(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        testSpace.globalDofs(element, testGdofs[elementIndex]);
        trialSpace.globalDofs(element, trialGdofs[elementIndex]);
        it->next();
    }

    // Distribute local matrices into the global matrix
    for (int e = 0; e < elementCount; ++e)
        for (int trialIndex = 0; trialIndex < trialGdofs[e].size(); ++trialIndex)
            for (int testIndex = 0; testIndex < testGdofs[e].size(); ++testIndex)
                result(testGdofs[e][testIndex], trialGdofs[e][trialIndex]) +=
                        localResult[e](testIndex, trialIndex);

    return std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >(
                new DiscreteDenseScalarValuedLinearOperator<ValueType>(result));
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
IdentityOperator<ValueType>::assembleWeakFormInSparseMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        typename IdentityOperator<ValueType>::LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
#ifdef WITH_TRILINOS
    // Fill local submatrices
    std::auto_ptr<GridView> view = testSpace.grid().leafView();
    const int elementCount = view->entityCount(0);
    std::vector<int> elementIndices(elementCount);
    for (int i = 0; i < elementCount; ++i)
        elementIndices[i] = i;
    std::vector<arma::Mat<ValueType> > localResult;
    assembler.evaluateLocalWeakForms(elementIndices, localResult);

    // Estimate number of entries in each row

//    This will be useful when we begin to use MPI
//    // Get global DOF indices for which this process is responsible
//    const int testGlobalDofCount = testSpace.globalDofCount();
//    Epetra_Map rowMap(testGlobalDofCount, 0 /* index-base */, comm);
//    std::vector<int> myTestGlobalDofs(rowMap.MyGlobalElements(),
//                                      rowMap.MyGlobalElements() +
//                                      rowMap.NumMyElements());
//    const int myTestGlobalDofCount = myTestGlobalDofs.size();

    const int testGlobalDofCount = testSpace.globalDofCount();
    arma::Col<int> nonzeroEntryCountEstimates(testGlobalDofCount);
    nonzeroEntryCountEstimates.fill(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > trialGdofs(elementCount);
    std::vector<std::vector<GlobalDofIndex> > testGdofs(elementCount);

    // Fill above lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        testSpace.globalDofs(element, testGdofs[elementIndex]);
        trialSpace.globalDofs(element, trialGdofs[elementIndex]);
        it->next();
    }

    // Upper estimate for the number of global trial DOFs coupled to a given
    // global test DOF: sum of the local trial DOF counts for each element that
    // contributes to the global test DOF in question
    for (int e = 0; e < elementCount; ++e)
        for (int testLdof = 0; testLdof < testGdofs[e].size(); ++testLdof)
            nonzeroEntryCountEstimates(testGdofs[e][testLdof]) +=
                    trialGdofs[e].size();

    Epetra_SerialComm comm; // To be replaced once we begin to use MPI
    Epetra_LocalMap rowMap(testGlobalDofCount, 0 /* index_base */, comm);
    std::auto_ptr<Epetra_FECrsMatrix> result(
                new Epetra_FECrsMatrix(Copy, rowMap,
                                       nonzeroEntryCountEstimates.memptr()));

    // TODO: make each process responsible for a subset of elements
    // Find maximum number of local dofs per element
    size_t maxLdofCount = 0;
    for (int e = 0; e < elementCount; ++e)
        maxLdofCount = std::max(maxLdofCount,
                                testGdofs[e].size() * trialGdofs[e].size());

    // Initialise sparse matrix with zeros at required positions
    arma::Col<double> zeros(maxLdofCount);
    zeros.fill(0.);
    for (int e = 0; e < elementCount; ++e)
        result->InsertGlobalValues(testGdofs[e].size(), &testGdofs[e][0],
                                  trialGdofs[e].size(), &trialGdofs[e][0],
                                  zeros.memptr());
    // Add contributions from individual elements
    for (int e = 0; e < elementCount; ++e)
        epetraSumIntoGlobalValues(
                    *result, testGdofs[e], trialGdofs[e], localResult[e]);
    result->GlobalAssemble();

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >(
                new DiscreteSparseScalarValuedLinearOperator<ValueType>(result));
#else // WITH_TRILINOS
    throw std::runtime_error("To enable assembly in sparse mode, recompile BEM++ "
                             "with the symbol WITH_TRILINOS defined.");
#endif
}

template <typename ValueType>
std::auto_ptr<typename IdentityOperator<ValueType>::LocalAssembler>
IdentityOperator<ValueType>::makeAssembler(
        const LocalAssemblerFactory& assemblerFactory,
        const GeometryFactory& geometryFactory,
        const Fiber::RawGridGeometry<ValueType>& rawGeometry,
        const std::vector<const Fiber::Basis<ValueType>*>& testBases,
        const std::vector<const Fiber::Basis<ValueType>*>& trialBases,
        const Fiber::OpenClHandler<ValueType, int>& openClHandler,
        bool /* cacheSingularIntegrals */) const
{
    return assemblerFactory.make(geometryFactory, rawGeometry,
                                 testBases, trialBases,
                                 m_expression, m_expression, this->multiplier(),
                                 openClHandler);
}

#ifdef COMPILE_FOR_FLOAT
template class IdentityOperator<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class IdentityOperator<double>;
#endif
// do we need this?
//#ifdef COMPILE_FOR_COMPLEX_FLOAT
//#include <complex>
//template class IdentityOperator<std::complex<float> >;
//#endif
//#ifdef COMPILE_FOR_COMPLEX_DOUBLE
//#include <complex>
//template class IdentityOperator<std::complex<double> >;
//#endif

} // namespace Bempp
