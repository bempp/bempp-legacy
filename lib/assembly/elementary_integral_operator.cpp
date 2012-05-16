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
#include "discrete_dense_linear_operator.hpp"
#include "evaluation_options.hpp"
#include "grid_function.hpp"
#include "interpolated_function.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/auto_timer.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/evaluator_for_integral_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include <armadillo>
#include <boost/make_shared.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <stdexcept>
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

namespace Bempp
{

// Body of parallel loop

namespace
{

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
    mutable typename Fiber::LocalAssemblerForOperators<ResultType>& m_assembler;
    // mutable OK because write access to this matrix is protected by a mutex
    mutable arma::Mat<ResultType>& m_result;

    // mutex must be mutable because we need to lock and unlock it
    mutable MutexType& m_mutex;
};

} // namespace

template <typename BasisFunctionType, typename KernelType, typename ResultType>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
ElementaryIntegralOperator(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace) :
    ElementaryLinearOperator<BasisFunctionType, ResultType>(testSpace, trialSpace)
{
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
bool
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (repr == AssemblyOptions::DENSE || repr == AssemblyOptions::ACA);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
makeAssembler(
        const LocalAssemblerFactory& assemblerFactory,
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
        const ParallelisationOptions& parallelisationOptions,
        bool cacheSingularIntegrals) const
{
    return assemblerFactory.makeAssemblerForIntegralOperators(
                testGeometryFactory, trialGeometryFactory,
                testRawGeometry, trialRawGeometry,
                testBases, trialBases,
                make_shared_from_ref(testExpression()),
                make_shared_from_ref(kernel()),
                make_shared_from_ref(trialExpression()),
                openClHandler, parallelisationOptions, cacheSingularIntegrals);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakForm(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    AutoTimer timer("\nAssembly took ");
    std::auto_ptr<LocalAssembler> assembler =
            makeAssemblerFromScratch(factory, options);
    return assembleWeakFormInternal(*assembler, options);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInternal(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    switch (options.operatorRepresentation()) {
    case AssemblyOptions::DENSE:
        return assembleWeakFormInDenseMode(assembler, options);
    case AssemblyOptions::ACA:
        return assembleWeakFormInAcaMode(assembler, options);
    case AssemblyOptions::FMM:
        throw std::runtime_error("ElementaryIntegralOperator::assembleWeakForm(): "
                                 "assembly mode FMM is not implemented yet");
    default:
        throw std::runtime_error("ElementaryIntegralOperator::assembleWeakForm(): "
                                 "invalid assembly mode");
    }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInDenseMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    const Space<BasisFunctionType>& testSpace = this->testSpace();
    const Space<BasisFunctionType>& trialSpace = this->trialSpace();

    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = trialSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > trialGlobalDofs(elementCount);
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        testSpace.globalDofs(element, testGlobalDofs[elementIndex]);
        trialSpace.globalDofs(element, trialGlobalDofs[elementIndex]);
        it->next();
    }

    // Make a vector of all element indices
    std::vector<int> testIndices(elementCount);
    for (int i = 0; i < elementCount; ++i)
        testIndices[i] = i;

    // Create the operator's matrix
    arma::Mat<ResultType> result(testSpace.globalDofCount(),
                                trialSpace.globalDofCount());
    result.fill(0.);

    typedef DenseWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
    typename Body::MutexType mutex;

    const ParallelisationOptions& parallelOptions =
            options.parallelisationOptions();
    int maxThreadCount = 1;
    if (parallelOptions.mode() == ParallelisationOptions::TBB) {
        if (parallelOptions.maxThreadCount() == ParallelisationOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = parallelOptions.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, elementCount),
                      Body(testIndices, testGlobalDofs, trialGlobalDofs,
                           assembler, result, mutex));

    //// Old serial code (TODO: decide whether to keep it behind e.g. #ifndef PARALLEL)
    //    std::vector<arma::Mat<ValueType> > localResult;
    //    // Loop over trial elements
    //    for (int trialIndex = 0; trialIndex < elementCount; ++trialIndex)
    //    {
    //        // Evaluate integrals over pairs of the current trial element and
    //        // all the test elements
    //        assembler.evaluateLocalWeakForms(TEST_TRIAL, testIndices, trialIndex,
    //                                         ALL_DOFS, localResult);

    //        // Loop over test indices
    //        for (int testIndex = 0; testIndex < elementCount; ++testIndex)
    //            // Add the integrals to appropriate entries in the operator's matrix
    //            for (int trialDof = 0; trialDof < trialGlobalDofs[trialIndex].size(); ++trialDof)
    //                for (int testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof)
    //                result(testGlobalDofs[testIndex][testDof],
    //                       trialGlobalDofs[trialIndex][trialDof]) +=
    //                        localResult[testIndex](testDof, trialDof);
    //    }

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteLinearOperator<ResultType> >(
                new DiscreteDenseLinearOperator<ResultType>(result));
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInAcaMode(
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    const Space<BasisFunctionType>& testSpace=this->testSpace();
    const Space<BasisFunctionType>& trialSpace=this->trialSpace();

    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleWeakForm(
                testSpace, trialSpace, assembler, options);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::Evaluator>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
makeEvaluator(
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const LocalAssemblerFactory& assemblerFactory,
        const EvaluationOptions& options) const
{
    // Collect the standard set of data necessary for construction of
    // evaluators and assemblers
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;
    typedef std::vector<std::vector<ResultType> > CoefficientsVector;
    typedef LocalAssemblerConstructionHelper Helper;

    shared_ptr<RawGridGeometry> rawGeometry;
    shared_ptr<GeometryFactory> geometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> bases;

    const Space<BasisFunctionType>& space = argument.space();
    Helper::collectGridData(space.grid(),
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelisationOptions().openClOptions(),
                              rawGeometry, openClHandler);
    Helper::collectBases(space, bases);

    // In addition, get coefficients of argument's expansion in each element
    const Grid& grid = space.grid();
    std::auto_ptr<GridView> view = grid.leafView();
    const int elementCount = view->entityCount(0);

    shared_ptr<CoefficientsVector> localCoefficients =
            boost::make_shared<CoefficientsVector>(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    for (int i = 0; i < elementCount; ++i) {
        const Entity<0>& element = it->entity();
        argument.getLocalCoefficients(element, (*localCoefficients)[i]);
        it->next();
    }

    // Now create the evaluator
    return assemblerFactory.makeEvaluatorForIntegralOperators(
                geometryFactory, rawGeometry,
                bases,
                make_shared_from_ref(kernel()),
                make_shared_from_ref(trialExpression()),
                localCoefficients,
                openClHandler);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<InterpolatedFunction<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
applyOffSurface(
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const Grid& evaluationGrid,
        const LocalAssemblerFactory& assemblerFactory,
        const EvaluationOptions& options) const
{
    std::auto_ptr<Evaluator> evaluator =
            makeEvaluator(argument, assemblerFactory, options);

    // Get coordinates of interpolation points, i.e. the evaluationGrid's vertices

    std::auto_ptr<GridView> evalView = evaluationGrid.leafView();
    const int evalGridDim = evaluationGrid.dim();
    const int evalPointCount = evalView->entityCount(evalGridDim);
    arma::Mat<CoordinateType> evalPoints(evalGridDim, evalPointCount);

    const IndexSet& evalIndexSet = evalView->indexSet();
    // TODO: extract into template function, perhaps add case evalGridDim == 1
    if (evalGridDim == 2) {
        const int vertexCodim = 2;
        std::auto_ptr<EntityIterator<vertexCodim> > it =
                evalView->entityIterator<vertexCodim>();
        while (!it->finished()) {
            const Entity<vertexCodim>& vertex = it->entity();
            const Geometry& geo = vertex.geometry();
            const int vertexIndex = evalIndexSet.entityIndex(vertex);
            arma::Col<CoordinateType> activeCol(evalPoints.unsafe_col(vertexIndex));
            geo.getCenter(activeCol);
            it->next();
        }
    } else if (evalGridDim == 3) {
        const int vertexCodim = 3;
        std::auto_ptr<EntityIterator<vertexCodim> > it =
                evalView->entityIterator<vertexCodim>();
        while (!it->finished()) {
            const Entity<vertexCodim>& vertex = it->entity();
            const Geometry& geo = vertex.geometry();
            const int vertexIndex = evalIndexSet.entityIndex(vertex);
            arma::Col<CoordinateType> activeCol(evalPoints.unsafe_col(vertexIndex));
            geo.getCenter(activeCol);
            it->next();
        }
    }

    // right now we don't bother about far and near field
    // (this might depend on evaluation options)

    arma::Mat<ResultType> result;
    evaluator->evaluate(Evaluator::FAR_FIELD, evalPoints, result);

    //    std::cout << "Interpolation results:\n";
    //    for (int point = 0; point < evalPointCount; ++point)
    //        std::cout << evalPoints(0, point) << "\t"
    //                  << evalPoints(1, point) << "\t"
    //                  << evalPoints(2, point) << "\t"
    //                  << result(0, point) << "\n";

    return std::auto_ptr<InterpolatedFunction<ResultType> >(
                new InterpolatedFunction<ResultType>(evaluationGrid, result));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryIntegralOperator);

} // namespace Bempp
