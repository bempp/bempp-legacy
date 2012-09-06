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

#include "evaluation_options.hpp"
#include "grid_function.hpp"
#include "interpolated_function.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../common/shared_ptr.hpp"

#include "../fiber/evaluator_for_integral_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/index_set.hpp"

namespace Bempp
{
template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<InterpolatedFunction<ResultType> >
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::evaluateOnGrid(
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const Grid& evaluationGrid,
        const QuadratureStrategy& quadStrategy,
        const EvaluationOptions& options) const
{
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
    std::auto_ptr<Evaluator> evaluator =
            makeEvaluator(argument, quadStrategy, options);

    // right now we don't bother about far and near field
    // (this might depend on evaluation options)

    arma::Mat<ResultType> result;
    evaluator->evaluate(Evaluator::FAR_FIELD, evaluationPoints, result);

    return result;
}

// UNDOCUMENTED PRIVATE METHODS

/** \cond PRIVATE */

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::Evaluator>
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
    Helper::collectGridData(space.grid(),
                            rawGeometry, geometryFactory);
    Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
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

/** \endcond */

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryPotentialOperator);

} // namespace Bempp
