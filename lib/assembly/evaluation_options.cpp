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

#include "evaluation_options.hpp"
#include "../common/global_parameters.hpp"

namespace Bempp {

EvaluationOptions::EvaluationOptions()
    : EvaluationOptions(GlobalParameters::parameterList()) {}

EvaluationOptions::EvaluationOptions(const ParameterList &parameters)
    : m_parameterList(parameters) {

  std::string assemblyType =
      parameters.get<std::string>("potentialOperatorAssemblyType");
  int maxThreadCount = parameters.get<int>("maxThreadCount");
  int verbosityLevel = parameters.get<int>("verbosityLevel");

  m_parallelizationOptions.setMaxThreadCount(maxThreadCount);

  if (assemblyType == "dense") {
    m_evaluationMode = DENSE;
  } else if (assemblyType == "hmat") {
    m_evaluationMode = HMAT;
  } else
    throw std::runtime_error(
        "EvaluationOptions::EvaluationOptions(): "
        "potentialoperatorAssemblyType has unsupported value.");

  if (verbosityLevel == -5)
    m_verbosityLevel = VerbosityLevel::LOW;
  else if (verbosityLevel == 0)
    m_verbosityLevel = VerbosityLevel::DEFAULT;
  else if (verbosityLevel == 5)
    m_verbosityLevel = VerbosityLevel::HIGH;
  else
    throw std::runtime_error(
        "Context::Context(): verbosityLevel has unsupported value");
}

EvaluationOptions::Mode EvaluationOptions::evaluationMode() const {
  return m_evaluationMode;
}

const AcaOptions &EvaluationOptions::acaOptions() const { return m_acaOptions; }

// void EvaluationOptions::switchToOpenCl(const OpenClOptions& openClOptions)
//{
//    m_parallelizationOptions.switchToOpenCl(openClOptions);
//}

void EvaluationOptions::setMaxThreadCount(int maxThreadCount) {
  m_parallelizationOptions.setMaxThreadCount(maxThreadCount);
}

void EvaluationOptions::switchToTbb(int maxThreadCount) {
  setMaxThreadCount(maxThreadCount);
}

const ParallelizationOptions &
EvaluationOptions::parallelizationOptions() const {
  return m_parallelizationOptions;
}

void EvaluationOptions::setVerbosityLevel(VerbosityLevel::Level level) {
  m_verbosityLevel = level;
}

VerbosityLevel::Level EvaluationOptions::verbosityLevel() const {
  return m_verbosityLevel;
}

} // namespace Bempp
