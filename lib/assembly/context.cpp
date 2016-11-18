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
Context<BasisFunctionType, ResultType>::Context()
    : Context(GlobalParameters::parameterList()) {}

template <typename BasisFunctionType, typename ResultType>
Context<BasisFunctionType, ResultType>::Context(
    const ParameterList &globalParameterList) {

  ParameterList parameters(globalParameterList);
  auto defaults = GlobalParameters::parameterList();

  std::string assemblyType = parameters.get<std::string>(
      "options.assembly.boundaryOperatorAssemblyType",
      defaults.get<std::string>(
          "options.assembly.boundaryOperatorAssemblyType"));
  if (assemblyType == "hmat")
    m_assemblyOptions.switchToHMatMode();
  else if (assemblyType == "dense")
    m_assemblyOptions.switchToDenseMode();
  else if (assemblyType == "cudadense")
    m_assemblyOptions.switchToCudaDenseMode();
  else
    throw std::runtime_error(
        "Context::Context(): boundaryOperatorAssemblyType has "
        "unsupported value.");

  m_assemblyOptions.setMaxThreadCount(
      parameters.get<int>("options.global.maxThreadCount",
                          defaults.get<int>("options.global.maxThreadCount")));

  int verbosityLevel =
      parameters.get<int>("options.global.verbosityLevel",
                          defaults.get<int>("options.global.verbosityLevel"));

  if (verbosityLevel == -5)
    m_assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
  else if (verbosityLevel == 0)
    m_assemblyOptions.setVerbosityLevel(VerbosityLevel::DEFAULT);
  else if (verbosityLevel == 5)
    m_assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH);
  else
    throw std::runtime_error(
        "Context::Context(): verbosityLevel has unsupported value");

  m_assemblyOptions.enableSingularIntegralCaching(parameters.get<bool>(
      "options.assembly.enableSingularIntegralCaching",
      defaults.get<bool>("options.assembly.enableSingularIntegralCaching")));

  m_assemblyOptions.enableBlasInQuadrature(AssemblyOptions::AUTO);

  // Cuda precision
  m_cudaOptions.setPrecision(
      parameters.get<std::string>("options.cuda.precision",
          defaults.get<std::string>("options.cuda.precision")));

  // Cuda element data caching
  bool elementDataCachingEnabled = parameters.get<bool>(
      "options.cuda.enableElementDataCaching",
      defaults.get<bool>("options.cuda.enableElementDataCaching"));
  if (elementDataCachingEnabled == true) {
    m_cudaOptions.enableElementDataCaching();
  } else {
    m_cudaOptions.disableElementDataCaching();
  }

  // Cuda kernel data caching
  bool kernelDataCachingEnabled = parameters.get<bool>(
      "options.cuda.enableKernelDataCaching",
      defaults.get<bool>("options.cuda.enableKernelDataCaching"));
  if (kernelDataCachingEnabled == true) {
    m_cudaOptions.enableKernelDataCaching();
  } else {
    m_cudaOptions.disableKernelDataCaching();
  }

  // Cuda device ids
//  m_cudaOptions.setDevices(
//      parameters.get<std::vector<int>>("options.cuda.deviceIds",
//          defaults.get<std::vector<int>>("options.cuda.deviceIds")));

  // Cuda numerical quadrature order
  m_cudaOptions.setQuadOrder(parameters.get<int>("options.cuda.quadOrder",
      defaults.get<int>("options.cuda.quadOrder")));

  // Cuda block size
  m_cudaOptions.setBlockSize(parameters.get<int>("options.cuda.blockSize",
      defaults.get<int>("options.cuda.blockSize")));

  // Cuda chunk size
  m_cudaOptions.setChunkElemPairCount(parameters.get<int>("options.cuda.chunkSize",
      defaults.get<int>("options.cuda.chunkSize")));

  Fiber::AccuracyOptionsEx accuracyOptions;

  accuracyOptions.setSingleRegular(
      parameters.get<double>(
          "options.quadrature.near.maxRelDist",
          defaults.get<double>("options.quadrature.near.maxRelDist")),
      parameters.get<int>(
          "options.quadrature.near.singleOrder",
          defaults.get<int>("options.quadrature.near.singleOrder")),
      parameters.get<double>(
          "options.quadrature.medium.maxRelDist",
          defaults.get<double>("options.quadrature.medium.maxRelDist")),
      parameters.get<int>(
          "options.quadrature.medium.singleOrder",
          defaults.get<int>("options.quadrature.medium.singleOrder")),
      parameters.get<int>(
          "options.quadrature.far.singleOrder",
          defaults.get<int>("options.quadrature.far.singleOrder")),
      false);

  accuracyOptions.setDoubleRegular(
      parameters.get<double>(
          "options.quadrature.near.maxRelDist",
          defaults.get<double>("options.quadrature.near.maxRelDist")),
      parameters.get<int>(
          "options.quadrature.near.doubleOrder",
          defaults.get<int>("options.quadrature.near.doubleOrder")),
      parameters.get<double>(
          "options.quadrature.medium.maxRelDist",
          defaults.get<double>("options.quadrature.medium.maxRelDist")),
      parameters.get<int>(
          "options.quadrature.medium.doubleOrder",
          defaults.get<int>("options.quadrature.medium.doubleOrder")),
      parameters.get<int>(
          "options.quadrature.far.doubleOrder",
          defaults.get<int>("options.quadrature.far.doubleOrder")),
      false);

  accuracyOptions.setDoubleSingular(
      parameters.get<int>(
          "options.quadrature.doubleSingular",
          defaults.get<int>("options.quadrature.doubleSingular")),
      false);

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
