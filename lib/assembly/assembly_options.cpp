#include "assembly_options.hpp"

#include <stdexcept>

namespace Bempp
{


AcaOptions::AcaOptions() :
    eps(1E-4),
    eta(0.8),
    maximumRank(10000),
    minimumBlockSize(16),
    recompress(true) {}



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

void AssemblyOptions::switchToSparse()
{
    m_representation = SPARSE;
}

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
