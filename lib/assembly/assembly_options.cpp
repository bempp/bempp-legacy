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

#include "assembly_options.hpp"

#include <stdexcept>

namespace Bempp {

AssemblyOptions::AssemblyOptions()
    : m_assemblyMode(DENSE), m_verbosityLevel(VerbosityLevel::DEFAULT),
      m_singularIntegralCaching(true), m_sparseStorageOfLocalOperators(true),
      m_jointAssembly(false), m_uniformQuadrature(true),
      m_blasInQuadrature(AUTO) {}

void AssemblyOptions::switchToDenseMode() { m_assemblyMode = DENSE; }

void AssemblyOptions::switchToHMatMode() { m_assemblyMode = HMAT; }

void AssemblyOptions::switchToAcaMode(const AcaOptions &acaOptions) {
  AcaOptions canonicalAcaOptions = acaOptions;
  if (!canonicalAcaOptions.globalAssemblyBeforeCompression) {
    canonicalAcaOptions.globalAssemblyBeforeCompression = true;
    canonicalAcaOptions.mode = AcaOptions::HYBRID_ASSEMBLY;
  }
  if ((int)canonicalAcaOptions.mode < AcaOptions::MIN_ASSEMBLY_MODE ||
      (int)canonicalAcaOptions.mode > AcaOptions::MAX_ASSEMBLY_MODE)
    throw std::invalid_argument("AssemblyOptions::switchToAcaMode(): "
                                "invalid ACA mode");
  if ((int)canonicalAcaOptions.reactionToUnsupportedMode <
          AcaOptions::MIN_REACTION ||
      (int)canonicalAcaOptions.reactionToUnsupportedMode >
          AcaOptions::MAX_REACTION)
    throw std::invalid_argument("AssemblyOptions::switchToAcaMode(): "
                                "invalid reaction to unsupported mode");
  m_assemblyMode = ACA;
  m_acaOptions = canonicalAcaOptions;
}

void AssemblyOptions::switchToDense() { switchToDenseMode(); }

void AssemblyOptions::switchToAca(const AcaOptions &acaOptions) {
  switchToAcaMode(acaOptions);
}

AssemblyOptions::Mode AssemblyOptions::assemblyMode() const {
  return m_assemblyMode;
}

const AcaOptions &AssemblyOptions::acaOptions() const { return m_acaOptions; }

// void AssemblyOptions::switchToOpenCl(const OpenClOptions& openClOptions)
//{
//    m_parallelizationOptions.switchToOpenCl(openClOptions);
//}

void AssemblyOptions::setMaxThreadCount(int maxThreadCount) {
  m_parallelizationOptions.setMaxThreadCount(maxThreadCount);
}

void AssemblyOptions::switchToTbb(int maxThreadCount) {
  setMaxThreadCount(maxThreadCount);
}

const ParallelizationOptions &AssemblyOptions::parallelizationOptions() const {
  return m_parallelizationOptions;
}

void AssemblyOptions::setVerbosityLevel(VerbosityLevel::Level level) {
  m_verbosityLevel = level;
}

VerbosityLevel::Level AssemblyOptions::verbosityLevel() const {
  return m_verbosityLevel;
}

void AssemblyOptions::enableSingularIntegralCaching(bool value) {
  m_singularIntegralCaching = value;
}

bool AssemblyOptions::isSingularIntegralCachingEnabled() const {
  return m_singularIntegralCaching;
}

void AssemblyOptions::enableSparseStorageOfLocalOperators(bool value) {
  m_sparseStorageOfLocalOperators = value;
}

bool AssemblyOptions::isSparseStorageOfLocalOperatorsEnabled() const {
  return m_sparseStorageOfLocalOperators;
}

void AssemblyOptions::enableSparseStorageOfMassMatrices(bool value) {
  enableSparseStorageOfLocalOperators(value);
}

bool AssemblyOptions::isSparseStorageOfMassMatricesEnabled() const {
  return isSparseStorageOfLocalOperatorsEnabled();
}

void AssemblyOptions::enableJointAssembly(bool value) {
  m_jointAssembly = value;
}

bool AssemblyOptions::isJointAssemblyEnabled() const { return m_jointAssembly; }

void AssemblyOptions::enableBlasInQuadrature(Value value) {
  if (value != AUTO && value != YES && value != NO)
    throw std::invalid_argument("AssemblyOptions::enableBlasInQuadrature(): "
                                "allowed values: AUTO, YES and NO");
  m_blasInQuadrature = value;
}

AssemblyOptions::Value AssemblyOptions::isBlasEnabledInQuadrature() const {
  return m_blasInQuadrature;
}

void AssemblyOptions::makeQuadratureOrderUniformInEachCluster(bool value) {
  m_uniformQuadrature = value;
}

bool AssemblyOptions::isQuadratureOrderUniformInEachCluster() const {
  return m_uniformQuadrature;
}

} // namespace Bempp
