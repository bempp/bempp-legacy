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
#include <limits>

namespace Bempp
{

AcaOptions::AcaOptions() :
    eps(1E-4),
    eta(1.2),
    minimumBlockSize(16),
    maximumBlockSize(std::numeric_limits<unsigned int>::max()),
    maximumRank(std::numeric_limits<unsigned int>::max()),
    recompress(true),
    outputPostscript(false),
    outputFname("aca.ps"),
    scaling(1.0)
{}



AssemblyOptions::AssemblyOptions() :
    m_representation(DENSE),
    m_singularIntegralCaching(AUTO)
{
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
    m_parallelisationOptions.switchToOpenCl(openClOptions);
}

void AssemblyOptions::switchToTbb(int maxThreadCount)
{
    m_parallelisationOptions.switchToTbb(maxThreadCount);
}

void AssemblyOptions::setSingularIntegralCaching(Mode mode)
{
    if (mode != AUTO && mode != NO && mode != YES)
        throw std::runtime_error("AssemblyOptions::setSingularIntegralCaching(): "
                                 "invalid mode");
    m_singularIntegralCaching = mode;
}

} // namespace Bempp
