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


#include "elementary_integral_operator.hpp"

#include "aca_global_assembler.hpp"
#include "assembly_options.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "context.hpp"
#include "evaluation_options.hpp"
#include "grid_function.hpp"
#include "interpolated_function.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/auto_timer.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/evaluator_for_integral_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/collection_of_basis_transformations.hpp"
#include "../fiber/quadrature_strategy.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/boost_ptr_vector_fwd.hpp"
#include <stdexcept>
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

namespace Bempp
{

// Helper functions and classes

namespace
{

// Body of parallel loop

template <typename BasisFunctionType, typename ResultType>
class DenseWeakFormAssemblerLoopBody
{
public:
    typedef tbb::spin_mutex MutexType;

    DenseWeakFormAssemblerLoopBody(
            const std::vector<int>& testIndices,
            const std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs,
            const std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs,
            Fiber::LocalAssemblerForOperators<ResultType>& assembler,
            arma::Mat<ResultType>& result, MutexType& mutex) :
        m_testIndices(testIndices),
        m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs),
        m_assembler(assembler), m_result(result), m_mutex(mutex) {
    }

    void operator() (const tbb::blocked_range<size_t>& r) const {
        const int elementCount = m_testIndices.size();
        std::vector<arma::Mat<ResultType> > localResult;
        for (size_t trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex) {
            // Evaluate integrals over pairs of the current trial element and
            // all the test elements
            m_assembler.evaluateLocalWeakForms(TEST_TRIAL, m_testIndices, trialIndex,
                                               ALL_DOFS, localResult);

            const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
            // Global assembly
            {
                MutexType::scoped_lock lock(m_mutex);
                // Loop over test indices
                for (int testIndex = 0; testIndex < elementCount; ++testIndex) {
                    const int testDofCount = m_testGlobalDofs[testIndex].size();
                    // Add the integrals to appropriate entries in the operator's matrix
                    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                        for (int testDof = 0; testDof < testDofCount; ++testDof)
                            m_result(m_testGlobalDofs[testIndex][testDof],
                                     m_trialGlobalDofs[trialIndex][trialDof]) +=
                                    localResult[testIndex](testDof, trialDof);
                }
            }
        }
    }

private:
    const std::vector<int>& m_testIndices;
    const std::vector<std::vector<GlobalDofIndex> >& m_testGlobalDofs;
    const std::vector<std::vector<GlobalDofIndex> >& m_trialGlobalDofs;
    // mutable OK because Assembler is thread-safe. (Alternative to "mutable" here:
    // make assembler's internal integrator map mutable)
    typename Fiber::LocalAssemblerForOperators<ResultType>& m_assembler;
    // mutable OK because write access to this matrix is protected by a mutex
    arma::Mat<ResultType>& m_result;

    // mutex must be mutable because we need to lock and unlock it
    MutexType& m_mutex;
};

/** Build a list of lists of global DOF indices corresponding to the local DOFs
 *  on each element of space.grid(). */
template <typename BasisFunctionType>
std::vector<std::vector<GlobalDofIndex> > gatherGlobalDofs(
        const Space<BasisFunctionType>& space)
{
    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = space.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > globalDofs(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        space.globalDofs(element, globalDofs[elementIndex]);
        it->next();
    }

    return globalDofs;
}

} // namespace

template <typename BasisFunctionType, typename KernelType, typename ResultType>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
ElementaryIntegralOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                           const shared_ptr<const Space<BasisFunctionType> >& range,
                           const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                           const std::string& label,
                           Symmetry symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
int
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
testComponentCount() const
{
    return testTransformations().argumentDimension();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
int
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
trialComponentCount() const
{
    return trialTransformations().argumentDimension();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
bool
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::isLocal() const
{
    return false;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
makeAssemblerImpl(
        const QuadratureStrategy& quadStrategy,
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        bool cacheSingularIntegrals) const
{
    return quadStrategy.makeAssemblerForIntegralOperators(
                testGeometryFactory, trialGeometryFactory,
                testRawGeometry, trialRawGeometry,
                testBases, trialBases,
                make_shared_from_ref(testTransformations()),
                make_shared_from_ref(kernels()),
                make_shared_from_ref(trialTransformations()),
                make_shared_from_ref(integral()),
                openClHandler, parallelizationOptions, cacheSingularIntegrals);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormImpl(
        const Context<BasisFunctionType, ResultType>& context) const
{
    AutoTimer timer("\nAssembly took ");
    std::auto_ptr<LocalAssembler> assembler =
            makeAssembler(context.quadStrategy(), context.assemblyOptions());
    return assembleWeakFormInternalImpl(*assembler, context.assemblyOptions());
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInternalImpl(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    switch (options.assemblyMode()) {
    case AssemblyOptions::DENSE:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInDenseMode(assembler, options).release());
    case AssemblyOptions::ACA:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInAcaMode(assembler, options).release());
    default:
        throw std::runtime_error(
                    "ElementaryIntegralOperator::assembleWeakFormInternalImpl(): "
                    "invalid assembly mode");
    }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInDenseMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs =
            gatherGlobalDofs(testSpace);
    std::vector<std::vector<GlobalDofIndex> > trialGlobalDofs =
            gatherGlobalDofs(trialSpace);
    const size_t testElementCount = testGlobalDofs.size();
    const size_t trialElementCount = trialGlobalDofs.size();

    // Make a vector of all element indices
    std::vector<int> testIndices(testElementCount);
    for (int i = 0; i < testElementCount; ++i)
        testIndices[i] = i;

    // Create the operator's matrix
    arma::Mat<ResultType> result(testSpace.globalDofCount(),
                                 trialSpace.globalDofCount());
    result.fill(0.);

    typedef DenseWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
    typename Body::MutexType mutex;

    const ParallelizationOptions& parallelOptions =
            options.parallelizationOptions();
    int maxThreadCount = 1;
    if (!parallelOptions.isOpenClEnabled()) {
        if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = parallelOptions.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount),
                      Body(testIndices, testGlobalDofs, trialGlobalDofs,
                           assembler, result, mutex));

    //// Old serial code (TODO: decide whether to keep it behind e.g. #ifndef PARALLEL)
    //    std::vector<arma::Mat<ValueType> > localResult;
    //    // Loop over trial elements
    //    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
    //    {
    //        // Evaluate integrals over pairs of the current trial element and
    //        // all the test elements
    //        assembler.evaluateLocalWeakForms(TEST_TRIAL, testIndices, trialIndex,
    //                                         ALL_DOFS, localResult);

    //        // Loop over test indices
    //        for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
    //            // Add the integrals to appropriate entries in the operator's matrix
    //            for (int trialDof = 0; trialDof < trialGlobalDofs[trialIndex].size(); ++trialDof)
    //                for (int testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof)
    //                result(testGlobalDofs[testIndex][testDof],
    //                       trialGlobalDofs[trialIndex][trialDof]) +=
    //                        localResult[testIndex](testDof, trialDof);
    //    }

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(
                new DiscreteDenseBoundaryOperator<ResultType>(result));
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInAcaMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
                testSpace, trialSpace, assembler, options,
                this->symmetry() & SYMMETRIC);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryIntegralOperator);

} // namespace Bempp
