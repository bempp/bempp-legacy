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

#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_trilinos.hpp"

#include "elementary_local_operator.hpp"

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#endif

#include "assembly_options.hpp"
#include "boundary_operator.hpp"
#include "cluster_construction_helper.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "context.hpp"

#include "../common/types.hpp"
#include "../common/complex_aux.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include <boost/type_traits/is_complex.hpp>

#include <tbb/tick_count.h>

#include <stdexcept>
#include <vector>

#ifdef WITH_TRILINOS
// This is a workaround of the problem of the abs() function being declared
// both in Epetra and in AHMED. It relies of the implementation detail (!) that
// in Epetra the declaration of abs is put between #ifndef __IBMCPP__ ...
// #endif. So it may well break in future versions of Trilinos. The ideal
// solution would be for AHMED to use namespaces.
#ifndef __IBMCPP__
#define __IBMCPP__
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>
#undef __IBMCPP__
#else
#include <Epetra_FECrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>
#endif
#endif // WITH_TRILINOS

namespace Bempp
{

namespace
{

#ifdef WITH_TRILINOS
// Internal helper functions for Epetra
template <typename ValueType>
int epetraSumIntoGlobalValues(Epetra_FECrsMatrix& matrix,
                              const std::vector<int>& rowIndices,
                              const std::vector<int>& colIndices,
                              const arma::Mat<ValueType>& values);

// Specialisation for double -- no intermediate array is needed
template <>
inline int epetraSumIntoGlobalValues<double>(
        Epetra_FECrsMatrix& matrix,
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

// Specialisation for float
template <>
inline int epetraSumIntoGlobalValues<float>(
        Epetra_FECrsMatrix& matrix,
        const std::vector<int>& rowIndices,
        const std::vector<int>& colIndices,
        const arma::Mat<float>& values)
{
    // Convert data from float into double (expected by Epetra)
    arma::Mat<double> doubleValues(values.n_rows, values.n_cols);
    std::copy(values.begin(), values.end(), doubleValues.begin());
    return epetraSumIntoGlobalValues<double>(
                matrix, rowIndices, colIndices, doubleValues);
}

// Specialisation for std::complex<float>.
// WARNING: at present only the real part is taken into account!
// This is sufficient as long as we provide real-valued basis functions only.
template <>
inline int epetraSumIntoGlobalValues<std::complex<float> >(
        Epetra_FECrsMatrix& matrix,
        const std::vector<int>& rowIndices,
        const std::vector<int>& colIndices,
        const arma::Mat<std::complex<float> >& values)
{
    // Extract the real part of "values" into an array of type double
    arma::Mat<double> doubleValues(values.n_rows, values.n_cols);
    // Check whether arma::real returns a view or a copy (if the latter,
    // this assert will fail)
    for (size_t i = 0; i < values.n_elem; ++i)
        doubleValues[i] = values[i].real();
    return epetraSumIntoGlobalValues<double>(
                matrix, rowIndices, colIndices, doubleValues);
}

// Specialisation for std::complex<double>.
// WARNING: at present only the real part is taken into account!
// This is sufficient as long as we provide real-valued basis functions only.
template <>
inline int epetraSumIntoGlobalValues<std::complex<double> >(
        Epetra_FECrsMatrix& matrix,
        const std::vector<int>& rowIndices,
        const std::vector<int>& colIndices,
        const arma::Mat<std::complex<double> >& values)
{
    // Extract the real part of "values" into an array of type double
    arma::Mat<double> doubleValues(values.n_rows, values.n_cols);
    for (size_t i = 0; i < values.n_elem; ++i)
        doubleValues[i] = values[i].real();
    return epetraSumIntoGlobalValues<double>(
                matrix, rowIndices, colIndices, doubleValues);
}
#endif

/** Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid(). */
template <typename BasisFunctionType>
void gatherGlobalDofs(
    const Space<BasisFunctionType>& testSpace,
    const Space<BasisFunctionType>& trialSpace,
    std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs,
    std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs,
    std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights,
    std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights)
{
    // We use the fact that test and trial space are required to be defined
    // on the same grid

    // Get the grid's leaf view so that we can iterate over elements
    const GridView& view = testSpace.gridView();
    const int elementCount = view.entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    testGlobalDofs.clear();
    testGlobalDofs.resize(elementCount);
    trialGlobalDofs.clear();
    trialGlobalDofs.resize(elementCount);
    // Weights of the local DOFs on elements
    testLocalDofWeights.clear();
    testLocalDofWeights.resize(elementCount);
    trialLocalDofWeights.clear();
    trialLocalDofWeights.resize(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view.elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        testSpace.getGlobalDofs(element, testGlobalDofs[elementIndex],
                                testLocalDofWeights[elementIndex]);
        trialSpace.getGlobalDofs(element, trialGlobalDofs[elementIndex],
                                 trialLocalDofWeights[elementIndex]);
        it->next();
    }
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
// ElementaryLocalOperator

template <typename BasisFunctionType, typename ResultType>
ElementaryLocalOperator<BasisFunctionType, ResultType>::ElementaryLocalOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
    if ((!domain->gridIsIdentical(*range)))
        throw std::invalid_argument(
                "ElementaryLocalOperator::ElementaryLocalOperator(): "
                "all three function spaces must be defined on the same grid.");
}

template <typename BasisFunctionType, typename ResultType>
bool ElementaryLocalOperator<BasisFunctionType, ResultType>::isLocal() const
{
    return true;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryLocalOperator<BasisFunctionType, ResultType>::assembleWeakFormImpl(
        const Context<BasisFunctionType, ResultType>& context) const
{
    bool verbose = (context.assemblyOptions().verbosityLevel() >=
                    VerbosityLevel::DEFAULT);
    if (verbose)
        std::cout << "Assembling the weak form of operator '"
                  << this->label() << "'..." << std::endl;

    tbb::tick_count start = tbb::tick_count::now();
    std::auto_ptr<LocalAssembler> assembler =this->makeAssembler(
                *context.quadStrategy(), context.assemblyOptions());
    shared_ptr<DiscreteBoundaryOperator<ResultType> > result =
            assembleWeakFormInternalImpl2(*assembler, context);
    tbb::tick_count end = tbb::tick_count::now();

    if (verbose)
        std::cout << "Assembly of the weak form of operator '" << this->label()
                  << "' took " << (end - start).seconds() << " s" << std::endl;
    return result;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryLocalOperator<BasisFunctionType, ResultType>::assembleWeakFormInternalImpl2(
        LocalAssembler& assembler,
        const Context<BasisFunctionType, ResultType>& context) const
{
#ifdef WITH_TRILINOS
    if (context.assemblyOptions().isSparseStorageOfMassMatricesEnabled())
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
            assembleWeakFormInSparseMode(assembler, context.assemblyOptions())
            .release());
#endif
    return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
        assembleWeakFormInDenseMode(assembler, context.assemblyOptions())
        .release());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryLocalOperator<BasisFunctionType, ResultType>::assembleWeakFormInDenseMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    // Fill local submatrices
    const GridView& view = testSpace.gridView();
    const size_t elementCount = view.entityCount(0);
    std::vector<int> elementIndices(elementCount);
    for (size_t i = 0; i < elementCount; ++i)
        elementIndices[i] = i;
    std::vector<arma::Mat<ResultType> > localResult;
    assembler.evaluateLocalWeakForms(elementIndices, localResult);

    // Create the operator's matrix
    arma::Mat<ResultType> result(testSpace.globalDofCount(),
                                 trialSpace.globalDofCount());
    result.fill(0.);

    // Retrieve global DOFs corresponding to local DOFs on all elements
    std::vector<std::vector<GlobalDofIndex> > testGdofs(elementCount);
    std::vector<std::vector<GlobalDofIndex> > trialGdofs(elementCount);
    std::vector<std::vector<BasisFunctionType> > testLdofWeights(elementCount);
    std::vector<std::vector<BasisFunctionType> > trialLdofWeights(elementCount);
    gatherGlobalDofs(testSpace, trialSpace, testGdofs, trialGdofs,
                     testLdofWeights, trialLdofWeights);

    // Distribute local matrices into the global matrix
    for (size_t e = 0; e < elementCount; ++e)
        for (size_t trialIndex = 0; trialIndex < trialGdofs[e].size(); ++trialIndex) {
            int trialGdof = trialGdofs[e][trialIndex];
            if (trialGdof < 0)
                continue;
            for (size_t testIndex = 0; testIndex < testGdofs[e].size(); ++testIndex) {
                int testGdof = testGdofs[e][testIndex];
                if (testGdof < 0)
                    continue;
                result(testGdof, trialGdof) +=
                        conj(testLdofWeights[e][testIndex]) *
                        trialLdofWeights[e][trialIndex] *
                        localResult[e](testIndex, trialIndex);
            }
        }

    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(
                new DiscreteDenseBoundaryOperator<ResultType>(result));
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryLocalOperator<BasisFunctionType, ResultType>::assembleWeakFormInSparseMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
#ifdef WITH_TRILINOS
    if (boost::is_complex<BasisFunctionType>::value)
        throw std::runtime_error(
                "ElementaryLocalOperator::assembleWeakFormInSparseMode(): "
                "sparse-mode assembly of identity operators for "
                "complex-valued basis functions is not supported yet");

    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    // Fill local submatrices
    const GridView& view = testSpace.gridView();
    const size_t elementCount = view.entityCount(0);
    std::vector<int> elementIndices(elementCount);
    for (size_t i = 0; i < elementCount; ++i)
        elementIndices[i] = i;
    std::vector<arma::Mat<ResultType> > localResult;
    assembler.evaluateLocalWeakForms(elementIndices, localResult);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGdofs(elementCount);
    std::vector<std::vector<GlobalDofIndex> > trialGdofs(elementCount);
    std::vector<std::vector<BasisFunctionType> > testLdofWeights(elementCount);
    std::vector<std::vector<BasisFunctionType> > trialLdofWeights(elementCount);
    gatherGlobalDofs(testSpace, trialSpace, testGdofs, trialGdofs,
                     testLdofWeights, trialLdofWeights);

    // Multiply matrix entries by DOF weights
    for (size_t e = 0; e < elementCount; ++e)
        for (size_t trialDof = 0; trialDof < trialGdofs[e].size(); ++trialDof)
            for (size_t testDof = 0; testDof < testGdofs[e].size(); ++testDof)
                localResult[e](testDof, trialDof) *=
                    conj(testLdofWeights[e][testDof]) *
                    trialLdofWeights[e][trialDof];

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
    const int trialGlobalDofCount = trialSpace.globalDofCount();
    arma::Col<int> nonzeroEntryCountEstimates(testGlobalDofCount);
    nonzeroEntryCountEstimates.fill(0);

    // Upper estimate for the number of global trial DOFs coupled to a given
    // global test DOF: sum of the local trial DOF counts for each element that
    // contributes to the global test DOF in question
    for (size_t e = 0; e < elementCount; ++e)
        for (size_t testLdof = 0; testLdof < testGdofs[e].size(); ++testLdof) {
            int testGdof = testGdofs[e][testLdof];
            if (testGdof >= 0)
                nonzeroEntryCountEstimates(testGdof) += trialGdofs[e].size();
        }

    Epetra_SerialComm comm; // To be replaced once we begin to use MPI
    Epetra_LocalMap rowMap(testGlobalDofCount, 0 /* index_base */, comm);
    Epetra_LocalMap colMap(trialGlobalDofCount, 0 /* index_base */, comm);
    shared_ptr<Epetra_FECrsMatrix> result = boost::make_shared<Epetra_FECrsMatrix>(
                Copy, rowMap, colMap,
                nonzeroEntryCountEstimates.memptr());

    // TODO: make each process responsible for a subset of elements
    // Find maximum number of local dofs per element
    size_t maxLdofCount = 0;
    for (size_t e = 0; e < elementCount; ++e)
        maxLdofCount = std::max(maxLdofCount,
                                testGdofs[e].size() * trialGdofs[e].size());

    // Initialise sparse matrix with zeros at required positions
    arma::Col<double> zeros(maxLdofCount);
    zeros.fill(0.);
    for (size_t e = 0; e < elementCount; ++e)
        result->InsertGlobalValues(testGdofs[e].size(), &testGdofs[e][0],
                                   trialGdofs[e].size(), &trialGdofs[e][0],
                                   zeros.memptr());
    // Add contributions from individual elements
    for (size_t e = 0; e < elementCount; ++e)
        epetraSumIntoGlobalValues(
            *result, testGdofs[e], trialGdofs[e], localResult[e]);
    result->GlobalAssemble();

    // If assembly mode is equal to ACA and we have AHMED,
    // construct the block cluster tree. Otherwise leave it uninitialized.
    typedef ClusterConstructionHelper<BasisFunctionType> CCH;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef ExtendedBemCluster<AhmedDofType> AhmedBemCluster;
    typedef bbxbemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;

    shared_ptr<AhmedBemBlcluster> blockCluster;
    shared_ptr<IndexPermutation> test_o2pPermutation, test_p2oPermutation;
    shared_ptr<IndexPermutation> trial_o2pPermutation, trial_p2oPermutation;
#ifdef WITH_AHMED
    if (options.assemblyMode() == AssemblyOptions::ACA) {
        const AcaOptions& acaOptions = options.acaOptions();
        bool indexWithGlobalDofs =
            acaOptions.mode != AcaOptions::HYBRID_ASSEMBLY;

        typedef ClusterConstructionHelper<BasisFunctionType> CCH;
        shared_ptr<AhmedBemCluster> testClusterTree;
        CCH::constructBemCluster(testSpace, indexWithGlobalDofs, acaOptions,
                                 testClusterTree,
                                 test_o2pPermutation, test_p2oPermutation);
        // TODO: construct a hermitian H-matrix if possible
        shared_ptr<AhmedBemCluster> trialClusterTree;
        CCH::constructBemCluster(trialSpace, indexWithGlobalDofs, acaOptions,
                                 trialClusterTree,
                                 trial_o2pPermutation, trial_p2oPermutation);
        unsigned int blockCount = 0;
        bool useStrongAdmissibilityCondition = !indexWithGlobalDofs;
        blockCluster.reset(CCH::constructBemBlockCluster(
                               acaOptions, false /* hermitian */,
                               *testClusterTree, *trialClusterTree,
                               useStrongAdmissibilityCondition,
                               blockCount)
                           .release());
    }
#endif

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(
        new DiscreteSparseBoundaryOperator<ResultType>(
                    result, this->symmetry(), NO_TRANSPOSE,
                    blockCluster, trial_o2pPermutation, test_o2pPermutation));
#else // WITH_TRILINOS
    throw std::runtime_error("ElementaryLocalOperator::assembleWeakFormInSparseMode(): "
                             "To enable assembly in sparse mode, recompile BEM++ "
                             "with the symbol WITH_TRILINOS defined.");
#endif
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<typename ElementaryLocalOperator<BasisFunctionType, ResultType>::LocalAssembler>
ElementaryLocalOperator<BasisFunctionType, ResultType>::makeAssemblerImpl(
        const QuadratureStrategy& quadStrategy,
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Fiber::Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<const Fiber::Shapeset<BasisFunctionType>*> >& trialShapesets,
        const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
        const ParallelizationOptions&,
        VerbosityLevel::Level /* verbosityLevel*/,
        bool /* cacheSingularIntegrals */) const
{
    if (testGeometryFactory.get() != trialGeometryFactory.get() ||
            testRawGeometry.get() != trialRawGeometry.get())
        throw std::invalid_argument("ElementaryLocalOperator::makeAssemblerImpl(): "
                                    "the test and trial spaces must be defined "
                                    "on the same grid");
    return quadStrategy.makeAssemblerForLocalOperators(
                testGeometryFactory, testRawGeometry,
                testShapesets, trialShapesets,
                make_shared_from_ref(testTransformations()),
                make_shared_from_ref(trialTransformations()),
                make_shared_from_ref(integral()),
                openClHandler);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ElementaryLocalOperator);

} // namespace Bempp
