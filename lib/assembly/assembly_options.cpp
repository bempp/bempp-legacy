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


#include "bempp/common/config_trilinos.hpp"
#include "assembly_options.hpp"

#include <stdexcept>
#include <limits>

namespace Bempp
{

AcaOptions::AcaOptions() :
    eps(1E-4),
    eta(1.2),
    minimumBlockSize(16),
    maximumBlockSize(std::numeric_limits<unsigned int>::max()),
    maximumRank(std::numeric_limits<unsigned int>::max()),
    globalAssemblyBeforeCompression(true),
    recompress(false),
    outputPostscript(false),
    outputFname("aca.ps"),
    scaling(1.0)
{
}



AssemblyOptions::AssemblyOptions() :
    m_assemblyMode(DENSE),
    m_verbosityLevel(VerbosityLevel::DEFAULT),
    m_singularIntegralCaching(true),
    m_sparseStorageOfMassMatrices(true)
{
}

void AssemblyOptions::switchToDenseMode()
{
    m_assemblyMode = DENSE;
}

void AssemblyOptions::switchToAcaMode(const AcaOptions& acaOptions)
{
    m_assemblyMode = ACA;
    m_acaOptions = acaOptions;
}

void AssemblyOptions::switchToDense()
{
    switchToDenseMode();
}

void AssemblyOptions::switchToAca(const AcaOptions& acaOptions)
{
    switchToAcaMode(acaOptions);
}

AssemblyOptions::Mode AssemblyOptions::assemblyMode() const {
    return m_assemblyMode;
}

const AcaOptions& AssemblyOptions::acaOptions() const {
    return m_acaOptions;
}

//void AssemblyOptions::switchToOpenCl(const OpenClOptions& openClOptions)
//{
//    m_parallelizationOptions.switchToOpenCl(openClOptions);
//}

void AssemblyOptions::setMaxThreadCount(int maxThreadCount)
{
    m_parallelizationOptions.setMaxThreadCount(maxThreadCount);
}

void AssemblyOptions::switchToTbb(int maxThreadCount)
{
    setMaxThreadCount(maxThreadCount);
}

const ParallelizationOptions& AssemblyOptions::parallelizationOptions() const
{
       return m_parallelizationOptions;
}

void AssemblyOptions::setVerbosityLevel(VerbosityLevel::Level level)
{
    m_verbosityLevel = level;
}

VerbosityLevel::Level AssemblyOptions::verbosityLevel() const
{
    return m_verbosityLevel;
}

void AssemblyOptions::enableSingularIntegralCaching(bool value)
{
    m_singularIntegralCaching = value;
}

bool AssemblyOptions::isSingularIntegralCachingEnabled() const
{
    return m_singularIntegralCaching;
}

void AssemblyOptions::enableSparseStorageOfMassMatrices(bool value)
{
    m_sparseStorageOfMassMatrices = value;
}

bool AssemblyOptions::isSparseStorageOfMassMatricesEnabled() const
{
    return m_sparseStorageOfMassMatrices;
}

} // namespace Bempp
