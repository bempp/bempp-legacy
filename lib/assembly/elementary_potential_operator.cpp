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

#include "hmat_global_assembler.hpp"
#include "assembled_potential_operator.hpp"
#include "evaluation_options.hpp"
#include "local_assembler_construction_helper.hpp"
#include "dense_global_assembler.hpp"
#include "discrete_boundary_operator.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"

#include "../fiber/evaluator_for_integral_operators.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/kernel_trial_integral.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"

#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/index_set.hpp"
#include "../grid/mapper.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
int ElementaryPotentialOperator<BasisFunctionType, KernelType,
                                ResultType>::componentCount() const {
  return integral().resultDimension();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
AssembledPotentialOperator<BasisFunctionType, ResultType>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
    assemble(const shared_ptr<const Space<BasisFunctionType>> &space,
             const shared_ptr<const Matrix<CoordinateType>> &evaluationPoints,
             const ParameterList &parameterList) const {
  if (!space)
    throw std::invalid_argument("ElementaryPotentialOperator::assemble(): "
                                "the shared pointer 'space' must not be null");
  if (!evaluationPoints)
    throw std::invalid_argument(
        "ElementaryPotentialOperator::assemble(): "
        "the shared pointer 'evaluationPoints' must not be null");
  if (evaluationPoints->rows() != space->grid()->dimWorld())
    throw std::invalid_argument(
        "ElementaryPotentialOperator::assemble(): "
        "the number of coordinates of each evaluation point must be "
        "equal to the dimension of the space containing the surface "
        "on which the function space 'space' is defined");

  auto quadStrategy =
      Context<BasisFunctionType, ResultType>(parameterList).quadStrategy();

  EvaluationOptions options(parameterList);

  std::unique_ptr<LocalAssembler> assembler =
      makeAssembler(*space, *evaluationPoints, *quadStrategy, options);
  shared_ptr<DiscreteBoundaryOperator<ResultType>> discreteOperator =
      assembleOperator(*space, *evaluationPoints /*TODO*/, *assembler,
                       parameterList);
  return AssembledPotentialOperator<BasisFunctionType, ResultType>(
      space, evaluationPoints, discreteOperator, componentCount());
}

// UNDOCUMENTED PRIVATE METHODS

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::unique_ptr<typename ElementaryPotentialOperator<
    BasisFunctionType, KernelType, ResultType>::LocalAssembler>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
    makeAssembler(const Space<BasisFunctionType> &space,
                  const Matrix<CoordinateType> &evaluationPoints,
                  const QuadratureStrategy &quadStrategy,
                  const EvaluationOptions &options) const {
  // Collect the standard set of data necessary for construction of
  // assemblers
  typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
  typedef std::vector<const Fiber::Shapeset<BasisFunctionType> *>
      ShapesetPtrVector;
  typedef LocalAssemblerConstructionHelper Helper;

  shared_ptr<RawGridGeometry> rawGeometry;
  shared_ptr<GeometryFactory> geometryFactory;
  shared_ptr<Fiber::OpenClHandler> openClHandler;
  shared_ptr<ShapesetPtrVector> shapesets;

  shared_ptr<const Grid> grid = space.grid();
  Helper::collectGridData(space, rawGeometry, geometryFactory);
  Helper::makeOpenClHandler(options.parallelizationOptions().openClOptions(),
                            rawGeometry, openClHandler);
  Helper::collectShapesets(space, shapesets);

  // Now create the assembler
  return quadStrategy.makeAssemblerForPotentialOperators(
      evaluationPoints, geometryFactory, rawGeometry, shapesets,
      make_shared_from_ref(kernels()),
      make_shared_from_ref(trialTransformations()),
      make_shared_from_ref(integral()), openClHandler,
      options.parallelizationOptions(), options.verbosityLevel());
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
    assembleOperator(const Space<BasisFunctionType> &space,
                     const Matrix<CoordinateType> &evaluationPoints,
                     LocalAssembler &assembler,
                     const ParameterList &parameterList) const {

  EvaluationOptions options(parameterList);

  switch (options.evaluationMode()) {
  case EvaluationOptions::DENSE:
    return shared_ptr<DiscreteBoundaryOperator<ResultType>>(
        assembleOperatorInDenseMode(space, evaluationPoints, assembler, options)
            .release());
  case EvaluationOptions::HMAT:
    return shared_ptr<DiscreteBoundaryOperator<ResultType>>(
        assembleOperatorInHMatMode(space, evaluationPoints, assembler,
                                   parameterList)
            .release());
  default:
    throw std::runtime_error(
        "ElementaryPotentialOperator::assembleWeakFormInternalImpl(): "
        "invalid assembly mode");
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
    assembleOperatorInDenseMode(const Space<BasisFunctionType> &space,
                                const Matrix<CoordinateType> &evaluationPoints,
                                LocalAssembler &assembler,
                                const EvaluationOptions &options) const {

  return DenseGlobalAssembler<BasisFunctionType, ResultType>::
      assemblePotentialOperator(evaluationPoints, space, assembler, options);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::unique_ptr<DiscreteBoundaryOperator<ResultType>>
ElementaryPotentialOperator<BasisFunctionType, KernelType, ResultType>::
    assembleOperatorInHMatMode(const Space<BasisFunctionType> &space,
                               const Matrix<CoordinateType> &evaluationPoints,
                               LocalAssembler &assembler,
                               const ParameterList &parameterList) const {
  return HMatGlobalAssembler<BasisFunctionType, ResultType>::
      assemblePotentialOperator(evaluationPoints, space, assembler,
                                parameterList);
}

/** \endcond */

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ElementaryPotentialOperator);

} // namespace Bempp
