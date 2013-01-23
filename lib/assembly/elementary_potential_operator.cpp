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

#include "elementary_potential_operator.hpp"

#include "aca_global_assembler.hpp"
#include "evaluation_options.hpp"
#include "grid_function.hpp"
#include "interpolated_function.hpp"
#include "local_assembler_construction_helper.hpp"
#include "discrete_null_boundary_operator.hpp"
#include "precalculated_potential_operator.hpp"

#include "../common/shared_ptr.hpp"

#include "../fiber/evaluator_for_integral_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/kernel_trial_integral.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/index_set.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
int
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
componentCount() const
{
    return integral().resultDimension();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<InterpolatedFunction<ResultType> >
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::evaluateOnGrid(
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const Grid& evaluationGrid,
        const QuadratureStrategy& quadStrategy,
        const EvaluationOptions& options) const
{
    if (evaluationGrid.dimWorld() != argument.grid()->dimWorld())
        throw std::invalid_argument(
                "ElementaryPotentialOperator::evaluateOnGrid(): "
                "the evaluation grid and the surface on which the grid "
                "function 'argument' is defined must be embedded in a space "
                "of the same dimension");

    std::auto_ptr<Evaluator> evaluator =
            makeEvaluator(argument, quadStrategy, options);

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

template <typename BasisFunctionType, typename KernelType, typename ResultType>
arma::Mat<ResultType>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::evaluateAtPoints(
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const arma::Mat<CoordinateType>& evaluationPoints,
        const QuadratureStrategy& quadStrategy,
        const EvaluationOptions& options) const
{
    if (evaluationPoints.n_rows != argument.grid()->dimWorld())
        throw std::invalid_argument(
                "ElementaryPotentialOperator::evaluateAtPoints(): "
                "the number of coordinates of each evaluation point must be "
                "equal to the dimension of the space containing the surface "
                "on which the grid function 'argument' is defined");

    std::auto_ptr<Evaluator> evaluator =
            makeEvaluator(argument, quadStrategy, options);

    // right now we don't bother about far and near field
    // (this might depend on evaluation options)

    arma::Mat<ResultType> result;
    evaluator->evaluate(Evaluator::FAR_FIELD, evaluationPoints, result);

    return result;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
//PrecalculatedPotentialOperator<BasisFunctionType, ResultType>
void
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
precalculateOperatorAtPoints(
        const Space<BasisFunctionType>& space,
        const arma::Mat<CoordinateType>& evaluationPoints,
        const QuadratureStrategy& quadStrategy,
        const AssemblyOptions& options)
{
    if (evaluationPoints.n_rows != space.grid()->dimWorld())
        throw std::invalid_argument(
                "ElementaryPotentialOperator::precalculateOperatorAtPoints(): "
                "the number of coordinates of each evaluation point must be "
                "equal to the dimension of the space containing the surface "
                "on which the function space 'space' is defined");

    std::auto_ptr<LocalAssembler> assembler =
            makeAssembler(space, evaluationPoints, quadStrategy, options);
    m_precalculatedOperator = assembleOperator(space, evaluationPoints,
                                               *assembler, options);

    //    shared_ptr<DiscreteBoundaryOperator<ResultType> > dop(
//                new DiscreteNullBoundaryOperator<ResultType>(5, 5));
//    return PrecalculatedPotentialOperator<BasisFunctionType, ResultType>(dop);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
arma::Mat<ResultType>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
applyPrecalculatedOperator(
        const GridFunction<BasisFunctionType, ResultType>& argument) const
{
    arma::Mat<ResultType> result(m_precalculatedOperator->rowCount(), 1);
    arma::Col<ResultType> colResult = result.unsafe_col(0);
    const arma::Col<ResultType>& coeffs = argument.coefficients();
    //std::cout << "Precalculated op\n" << m_precalculatedOperator->asMatrix();
    m_precalculatedOperator->apply(NO_TRANSPOSE, coeffs, colResult, 1., 0.);
    //std::cout << "colResult\n" << colResult;
    int componentCount = integral().resultDimension();
    assert(result.n_rows % componentCount == 0);
    result.reshape(componentCount, result.n_rows / componentCount);
    return result;
}

// UNDOCUMENTED PRIVATE METHODS

/** \cond PRIVATE */

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryPotentialOperator<
BasisFunctionType, KernelType, ResultType>::Evaluator>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::makeEvaluator(
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const QuadratureStrategy& quadStrategy,
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

    const Space<BasisFunctionType>& space = *argument.space();
    shared_ptr<const Grid> grid = space.grid();
    Helper::collectGridData(*grid,
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                              rawGeometry, openClHandler);
    Helper::collectBases(space, bases);

    // In addition, get coefficients of argument's expansion in each element
    std::auto_ptr<GridView> view = grid->leafView();
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
    return quadStrategy.makeEvaluatorForIntegralOperators(
                geometryFactory, rawGeometry,
                bases,
                make_shared_from_ref(kernels()),
                make_shared_from_ref(trialTransformations()),
                make_shared_from_ref(integral()),
                localCoefficients,
                openClHandler,
                options.parallelizationOptions());
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryPotentialOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::makeAssembler(
        const Space<BasisFunctionType>& space,
        const arma::Mat<CoordinateType>& evaluationPoints,
        const QuadratureStrategy& quadStrategy,
        const AssemblyOptions& options) const
{
    // Collect the standard set of data necessary for construction of
    // assemblers
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;
    typedef std::vector<std::vector<ResultType> > CoefficientsVector;
    typedef LocalAssemblerConstructionHelper Helper;

    shared_ptr<RawGridGeometry> rawGeometry;
    shared_ptr<GeometryFactory> geometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> bases;

    shared_ptr<const Grid> grid = space.grid();
    Helper::collectGridData(*grid,
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                              rawGeometry, openClHandler);
    Helper::collectBases(space, bases);

    // Now create the evaluator
    return quadStrategy.makeAssemblerForPotentialOperators(
                evaluationPoints,
                geometryFactory, rawGeometry,
                bases,
                make_shared_from_ref(kernels()),
                make_shared_from_ref(trialTransformations()),
                make_shared_from_ref(integral()),
                openClHandler,
                options.parallelizationOptions(),
                options.verbosityLevel());
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
assembleOperator(
        const Space<BasisFunctionType>& space,
        const arma::Mat<CoordinateType>& evaluationPoints,
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    switch (options.assemblyMode()) {
    case AssemblyOptions::DENSE:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleOperatorInDenseMode(space, evaluationPoints,
                                                assembler, options).release());
    case AssemblyOptions::ACA:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleOperatorInAcaMode(space, evaluationPoints,
                                              assembler, options).release());
    default:
        throw std::runtime_error(
                    "ElementaryPotentialOperator::assembleWeakFormInternalImpl(): "
                    "invalid assembly mode");
    }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
assembleOperatorInDenseMode(
        const Space<BasisFunctionType>& space,
        const arma::Mat<CoordinateType>& evaluationPoints,
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >();
//    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
//    const Space<BasisFunctionType>& trialSpace = *this->domain();

//    // Global DOF indices corresponding to local DOFs on elements
//    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs, trialGlobalDofs;
//    std::vector<std::vector<BasisFunctionType> > testLocalDofWeights,
//        trialLocalDofWeights;
//    gatherGlobalDofs(testSpace, testGlobalDofs, testLocalDofWeights);
//    if (&testSpace == &trialSpace) {
//        trialGlobalDofs = testGlobalDofs;
//        trialLocalDofWeights = testLocalDofWeights;
//    } else
//        gatherGlobalDofs(trialSpace, trialGlobalDofs, trialLocalDofWeights);
//    const size_t testElementCount = testGlobalDofs.size();
//    const size_t trialElementCount = trialGlobalDofs.size();

//    // Make a vector of all element indices
//    std::vector<int> testIndices(testElementCount);
//    for (int i = 0; i < testElementCount; ++i)
//        testIndices[i] = i;

//    // Create the operator's matrix
//    arma::Mat<ResultType> result(testSpace.globalDofCount(),
//                                 trialSpace.globalDofCount());
//    result.fill(0.);

//    typedef DenseWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
//    typename Body::MutexType mutex;

//    const ParallelizationOptions& parallelOptions =
//            options.parallelizationOptions();
//    int maxThreadCount = 1;
//    if (!parallelOptions.isOpenClEnabled()) {
//        if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
//            maxThreadCount = tbb::task_scheduler_init::automatic;
//        else
//            maxThreadCount = parallelOptions.maxThreadCount();
//    }
//    tbb::task_scheduler_init scheduler(maxThreadCount);
//    {
//        Fiber::SerialBlasRegion region;
//        tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount),
//                          Body(testIndices, testGlobalDofs, trialGlobalDofs,
//                               testLocalDofWeights, trialLocalDofWeights,
//                               assembler, result, mutex));
//    }

//    //// Old serial code (TODO: decide whether to keep it behind e.g. #ifndef PARALLEL)
//    //    std::vector<arma::Mat<ValueType> > localResult;
//    //    // Loop over trial elements
//    //    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
//    //    {
//    //        // Evaluate integrals over pairs of the current trial element and
//    //        // all the test elements
//    //        assembler.evaluateLocalWeakForms(TEST_TRIAL, testIndices, trialIndex,
//    //                                         ALL_DOFS, localResult);

//    //        // Loop over test indices
//    //        for (int testIndex = 0; testIndex < testElementCount; ++testIndex)
//    //            // Add the integrals to appropriate entries in the operator's matrix
//    //            for (int trialDof = 0; trialDof < trialGlobalDofs[trialIndex].size(); ++trialDof)
//    //                for (int testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof)
//    //                result(testGlobalDofs[testIndex][testDof],
//    //                       trialGlobalDofs[trialIndex][trialDof]) +=
//    //                        localResult[testIndex](testDof, trialDof);
//    //    }

//    // Create and return a discrete operator represented by the matrix that
//    // has just been calculated
//    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(
//                new DiscreteDenseBoundaryOperator<ResultType>(result));
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
assembleOperatorInAcaMode(
        const Space<BasisFunctionType>& space,
        const arma::Mat<CoordinateType>& evaluationPoints,
        LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
                evaluationPoints, space, assembler, options);
}

/** \endcond */

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryPotentialOperator);

} // namespace Bempp
