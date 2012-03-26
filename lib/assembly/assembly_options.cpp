#include "assembly_options.hpp"

#include <stdexcept>

namespace Bempp
{

AssemblyOptions::AssemblyOptions() :
    // TODO: perhaps set m_acaOptions to some defaults
    m_representation(DENSE),
    m_parallelism(TBB), m_maxThreadCount(AUTO),
    m_singularIntegralCaching(AUTO)
{}

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
}

void AssemblyOptions::switchToTbb(int maxThreadCount)
{
    if (maxThreadCount <= 0 && maxThreadCount != AUTO)
        throw std::runtime_error("AssemblyOptions::switchToTbb(): "
                                 "maxThreadCount must be positive or equal to AUTO");
    m_parallelism = TBB;
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
