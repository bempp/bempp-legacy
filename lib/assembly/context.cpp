// Copyright (C) 2011 by the BEM++ Authors
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

#include "context.hpp"

#include "abstract_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/verbosity_level.hpp"
#include "../fiber/accuracy_options.hpp"
#include "numerical_quadrature_strategy.hpp"
#include <Teuchos_ParameterList.hpp>

#include <boost/make_shared.hpp>
#include <iostream>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context(
    const shared_ptr<const QuadratureStrategy> &quadStrategy,
    const AssemblyOptions &assemblyOptions,
    const ParameterList &globalParameterList)
    : m_quadStrategy(quadStrategy), m_assemblyOptions(assemblyOptions),
      m_globalParameterList(globalParameterList) {
  if (quadStrategy.get() == 0)
    throw std::invalid_argument("Context::Context(): "
                                "quadStrategy must not be null");
}

template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context() :
    Context(GlobalParameters::parameterList()){}


template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context(
    const ParameterList &globalParameterList) {

  ParameterList parameters(globalParameterList);
  parameters.setParametersNotAlreadySet(GlobalParameters::parameterList());

  std::string assemblyType =
      parameters.get<std::string>("boundaryOperatorAssemblyType");
  if (assemblyType == "hmat")
    m_assemblyOptions.switchToHMatMode();
  else if (assemblyType == "dense")
    m_assemblyOptions.switchToDenseMode();
  else
    throw std::runtime_error(
        "Context::Context(): boundaryOperatorAssemblyType has "
        "unsupported value.");

  m_assemblyOptions.setMaxThreadCount(parameters.get<int>("maxThreadCount"));

  int verbosityLevel = parameters.get<int>("verbosityLevel");

  if (verbosityLevel == -5)
    m_assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
  else if (verbosityLevel == 0)
    m_assemblyOptions.setVerbosityLevel(VerbosityLevel::DEFAULT);
  else if (verbosityLevel == 5)
    m_assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH);
  else
    throw std::runtime_error(
        "Context::Context(): verbosityLevel has unsupported value");

  m_assemblyOptions.enableSingularIntegralCaching(
      parameters.get<bool>("enableSingularIntegralCaching"));

  std::string enableBlasInQuadrature =
      parameters.get<std::string>("enableBlasInQuadrature");
  if (enableBlasInQuadrature == "auto")
    m_assemblyOptions.enableBlasInQuadrature(AssemblyOptions::AUTO);
  else if (enableBlasInQuadrature == "yes")
    m_assemblyOptions.enableBlasInQuadrature(AssemblyOptions::YES);
  else if (enableBlasInQuadrature == "no")
    m_assemblyOptions.enableBlasInQuadrature(AssemblyOptions::NO);
  else
    throw std::runtime_error("Context::Context(): enableBlasInQuadrature "
                             "has unsupported value");

  Fiber::AccuracyOptionsEx accuracyOptions;

  auto quadOps = parameters.sublist("QuadratureOrders");

  accuracyOptions.setSingleRegular(
      quadOps.sublist("near").get<double>("maxRelDist"),
      quadOps.sublist("near").get<int>("singleOrder"),
      quadOps.sublist("medium").get<double>("maxRelDist"),
      quadOps.sublist("medium").get<int>("singleOrder"),
      quadOps.sublist("far").get<int>("singleOrder"),
      quadOps.get<bool>("quadratureOrdersAreRelative"));

  accuracyOptions.setDoubleRegular(
      quadOps.sublist("near").get<double>("maxRelDist"),
      quadOps.sublist("near").get<int>("doubleOrder"),
      quadOps.sublist("medium").get<double>("maxRelDist"),
      quadOps.sublist("medium").get<int>("doubleOrder"),
      quadOps.sublist("far").get<int>("doubleOrder"),
      quadOps.get<bool>("quadratureOrdersAreRelative"));

  accuracyOptions.setDoubleSingular(
      quadOps.get<int>("doubleSingular"),
      quadOps.get<bool>("quadratureOrdersAreRelative"));

  m_quadStrategy.reset(
      new NumericalQuadratureStrategy<BasisFunctionType, ResultType>(
          accuracyOptions));

  m_globalParameterList = parameters;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
Context<BasisFunctionType, ResultType>::getWeakForm(
    const AbstractBoundaryOperator<BasisFunctionType, ResultType> &op) const {
  return op.assembleWeakForm(*this);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

} // namespace Bempp
