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


#include "config_trilinos.hpp"
#include "assembly_options.hpp"

#include <stdexcept>

namespace Bempp
{


AcaOptions::AcaOptions() :
    eps(1E-4),
    eta(1.2),
    minimumBlockSize(16),
    maximumRank(10000),
    recompress(true),
    outputPostscript(false),
    scaling(1.0),
    outputFname("aca.ps")
{}



AssemblyOptions::AssemblyOptions() :
    // TODO: perhaps set m_acaOptions to some defaults
    m_representation(DENSE),
    m_parallelism(TBB), m_maxThreadCount(AUTO),
    m_singularIntegralCaching(AUTO)
{
    m_openClOptions.useOpenCl = false;
}

void AssemblyOptions::switchToDense()
{
    m_representation = DENSE;
}

void AssemblyOptions::switchToAca(const AcaOptions& acaOptions)
{
    m_representation = ACA;
    m_acaOptions = acaOptions;
}

void AssemblyOptions::switchToFmm()
{
    m_representation = FMM;
}

//void AssemblyOptions::switchToSparse()
//{
//    m_representation = SPARSE;
//}

void AssemblyOptions::switchToOpenCl(const OpenClOptions& openClOptions)
{
    m_parallelism = OPEN_CL;
    m_openClOptions = openClOptions;
    m_openClOptions.useOpenCl = true;
}

void AssemblyOptions::switchToTbb(int maxThreadCount)
{
    m_parallelism = TBB;
    m_openClOptions.useOpenCl = false;
    if (maxThreadCount <= 0 && maxThreadCount != AUTO)
        throw std::runtime_error("AssemblyOptions::switchToTbb(): "
                                 "maxThreadCount must be positive or equal to AUTO");
    m_maxThreadCount = maxThreadCount;
}

void AssemblyOptions::setSingularIntegralCaching(Mode mode)
{
    if (mode != AUTO && mode != NO && mode != YES)
        throw std::runtime_error("AssemblyOptions::setSingularIntegralCaching(): "
                                 "invalid mode");
    m_singularIntegralCaching = mode;
}

} // namespace Bempp
