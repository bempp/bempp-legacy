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


#include "bempp/common/config_trilinos.hpp"

#include "identity_operator.hpp"

#include "assembly_options.hpp"
#include "boundary_operator.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "context.hpp"

#include "../common/auto_timer.hpp"
#include "../common/types.hpp"
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
#include <stdexcept>
#include <vector>

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

} // anonymous namespace
#endif

////////////////////////////////////////////////////////////////////////////////
// IdentityOperatorId

template <typename BasisFunctionType, typename ResultType>
IdentityOperatorId<BasisFunctionType, ResultType>::IdentityOperatorId(
        const IdentityOperator<BasisFunctionType, ResultType>& op) :
    m_domain(op.domain().get()), m_range(op.range().get()),
    m_dualToRange(op.dualToRange().get())
{
}

template <typename BasisFunctionType, typename ResultType>
size_t IdentityOperatorId<BasisFunctionType, ResultType>::hash() const
{
    typedef IdentityOperator<BasisFunctionType, ResultType>
            OperatorType;
    size_t result = tbb::tbb_hasher(typeid(OperatorType).name());
    tbb_hash_combine(result, m_domain);
    tbb_hash_combine(result, m_range);
    tbb_hash_combine(result, m_dualToRange);
    return result;
}

template <typename BasisFunctionType, typename ResultType>
bool IdentityOperatorId<BasisFunctionType, ResultType>::isEqual(
        const AbstractBoundaryOperatorId &other) const
{
    // dynamic_cast won't suffice since we want to make sure both objects
    // are of exactly the same type (dynamic_cast would succeed for a subclass)
    if (typeid(other) == typeid(*this))
    {
        const IdentityOperatorId& otherCompatible =
            static_cast<const IdentityOperatorId&>(other);
        return (m_domain == otherCompatible.m_domain &&
                m_range == otherCompatible.m_range &&
                m_dualToRange == otherCompatible.m_dualToRange);
    }
    else
        return false;
}

////////////////////////////////////////////////////////////////////////////////
// IdentityOperator

template <typename BasisFunctionType, typename ResultType>
struct IdentityOperator<BasisFunctionType, ResultType>::Impl
{
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;

    Impl() : transformations(TransformationFunctor())
    {}

    Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
};

template <typename BasisFunctionType, typename ResultType>
IdentityOperator<BasisFunctionType, ResultType>::IdentityOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label) :
    Base(domain, range, dualToRange, label,
         domain == dualToRange ? HERMITIAN : NO_SYMMETRY),
    m_impl(new Impl),
    m_id(boost::make_shared<IdentityOperatorId<BasisFunctionType, ResultType> >(
             *this))
{
}

template <typename BasisFunctionType, typename ResultType>
IdentityOperator<BasisFunctionType, ResultType>::IdentityOperator(
        const IdentityOperator& other) :
    Base(other), m_impl(new Impl(*other.m_impl)),
    m_id(other.m_id)
{
}

template <typename BasisFunctionType, typename ResultType>
IdentityOperator<BasisFunctionType, ResultType>::~IdentityOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
IdentityOperator<BasisFunctionType, ResultType>::id() const
{
    return m_id;
}

template <typename BasisFunctionType, typename ResultType>
bool IdentityOperator<BasisFunctionType, ResultType>::isLocal() const
{
    return true;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
IdentityOperator<BasisFunctionType, ResultType>::assembleWeakFormImpl(
        const Context<BasisFunctionType, ResultType>& context) const
{
    AutoTimer timer("\nAssembly took ");
    std::auto_ptr<LocalAssembler> assembler = makeAssembler(
                context.quadStrategy(), context.assemblyOptions());
    return assembleWeakFormInternalImpl(*assembler, context.assemblyOptions());
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
IdentityOperator<BasisFunctionType, ResultType>::assembleWeakFormInternalImpl(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
#ifdef WITH_TRILINOS
    if (options.isSparseStorageOfMassMatricesEnabled())
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInSparseMode(assembler, options).release());
#endif
    return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInDenseMode(assembler, options).release());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
IdentityOperator<BasisFunctionType, ResultType>::assembleWeakFormInDenseMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    // Fill local submatrices
    std::auto_ptr<GridView> view = testSpace.grid().leafView();
    const size_t elementCount = view->entityCount(0);
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
    for (size_t e = 0; e < elementCount; ++e)
        for (size_t trialIndex = 0; trialIndex < trialGdofs[e].size(); ++trialIndex)
            for (size_t testIndex = 0; testIndex < testGdofs[e].size(); ++testIndex)
                result(testGdofs[e][testIndex], trialGdofs[e][trialIndex]) +=
                        localResult[e](testIndex, trialIndex);

    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(
                new DiscreteDenseBoundaryOperator<ResultType>(result));
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
IdentityOperator<BasisFunctionType, ResultType>::assembleWeakFormInSparseMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
#ifdef WITH_TRILINOS
    if (boost::is_complex<BasisFunctionType>::value)
        throw std::runtime_error(
                "IdentityOperator::assembleWeakFormInSparseMode(): "
                "sparse-mode assembly of identity operators for "
                "complex-valued basis functions is not supported yet");

    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    // Fill local submatrices
    std::auto_ptr<GridView> view = testSpace.grid().leafView();
    const size_t elementCount = view->entityCount(0);
    std::vector<int> elementIndices(elementCount);
    for (size_t i = 0; i < elementCount; ++i)
        elementIndices[i] = i;
    std::vector<arma::Mat<ResultType> > localResult;
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
    const int trialGlobalDofCount = trialSpace.globalDofCount();
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
    for (size_t e = 0; e < elementCount; ++e)
        for (size_t testLdof = 0; testLdof < testGdofs[e].size(); ++testLdof)
            nonzeroEntryCountEstimates(testGdofs[e][testLdof]) +=
                    trialGdofs[e].size();

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

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(
                new DiscreteSparseBoundaryOperator<ResultType>(result));
#else // WITH_TRILINOS
    throw std::runtime_error("IdentityOperator::assembleWeakFormInSparseMode(): "
                             "To enable assembly in sparse mode, recompile BEM++ "
                             "with the symbol WITH_TRILINOS defined.");
#endif
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<typename IdentityOperator<BasisFunctionType, ResultType>::LocalAssembler>
IdentityOperator<BasisFunctionType, ResultType>::makeAssemblerImpl(
        const QuadratureStrategy& quadStrategy,        
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
        const ParallelizationOptions&,
        bool /* cacheSingularIntegrals */) const
{
    shared_ptr<const Fiber::CollectionOfBasisTransformations<CoordinateType> >
            transformations = make_shared_from_ref(m_impl->transformations);

    if (testGeometryFactory.get() != trialGeometryFactory.get() ||
            testRawGeometry.get() != trialRawGeometry.get())
        throw std::invalid_argument("IdentityOperator::makeAssemblerImpl(): "
                                    "the test and trial spaces must be defined "
                                    "on the same grid");
    return quadStrategy.makeAssemblerForIdentityOperators(
                testGeometryFactory, testRawGeometry,
                testBases, trialBases,
                transformations, transformations,
                openClHandler);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
identityOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& domain,
                 const shared_ptr<const Space<BasisFunctionType> >& range,
                 const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                 const std::string& label)
{
    typedef IdentityOperator<BasisFunctionType, ResultType> Id;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                context,
                boost::make_shared<Id>(domain, range, dualToRange, label));
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    identityOperator( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const std::string&)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(IdentityOperator);

} // namespace Bempp
