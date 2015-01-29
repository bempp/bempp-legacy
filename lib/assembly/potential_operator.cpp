#include "potential_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "assembled_potential_operator.hpp"
#include "interpolated_function.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
AssembledPotentialOperator<BasisFunctionType, ResultType>
PotentialOperator<BasisFunctionType, ResultType>::assemble(
    const shared_ptr<const Space<BasisFunctionType>> &space,
    const shared_ptr<const arma::Mat<CoordinateType>> &evaluationPoints,
    const ParameterList &parameterList) const {

  shared_ptr<Context<BasisFunctionType, ResultType>> context(
      new Context<BasisFunctionType, ResultType>(parameterList));

  return this->assemble(space, evaluationPoints, *context->quadStrategy(),
                        EvaluationOptions(parameterList));
}

template <typename BasisFunctionType, typename ResultType>
arma::Mat<ResultType>
PotentialOperator<BasisFunctionType, ResultType>::evaluateAtPoints(
    const GridFunction<BasisFunctionType, ResultType> &argument,
    const arma::Mat<CoordinateType> &evaluationPoints,
    const ParameterList &parameterList) const {

  shared_ptr<Context<BasisFunctionType, ResultType>> context(
      new Context<BasisFunctionType, ResultType>(parameterList));

  return this->evaluateAtPoints(argument, evaluationPoints,
                                *context->quadStrategy(),
                                EvaluationOptions(parameterList));
}

template <typename BasisFunctionType, typename ResultType>
std::unique_ptr<InterpolatedFunction<ResultType>>
PotentialOperator<BasisFunctionType, ResultType>::evaluateOnGrid(
    const GridFunction<BasisFunctionType, ResultType> &argument,
    const Grid &evaluationGrid, const ParameterList &parameterList) const {

  shared_ptr<Context<BasisFunctionType, ResultType>> context(
      new Context<BasisFunctionType, ResultType>(parameterList));

  return this->evaluateOnGrid(argument, evaluationGrid,
                              *context->quadStrategy(),
                              EvaluationOptions(parameterList));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
        PotentialOperator);

}

