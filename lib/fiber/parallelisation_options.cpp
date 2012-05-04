#include "parallelisation_options.hpp"

#include <stdexcept>

namespace Fiber
{

ParallelisationOptions::ParallelisationOptions() :
    m_mode(TBB), m_maxThreadCount(AUTO)
{
    m_openClOptions.useOpenCl = false;
}

void ParallelisationOptions::switchToOpenCl(const OpenClOptions& openClOptions)
{
    m_mode = OPEN_CL;
    m_openClOptions = openClOptions;
    m_openClOptions.useOpenCl = true;
}

void ParallelisationOptions::switchToTbb(int maxThreadCount)
{
    m_mode = TBB;
    m_openClOptions.useOpenCl = false;
    if (maxThreadCount <= 0 && maxThreadCount != AUTO)
        throw std::runtime_error("ParallelisationOptions::switchToTbb(): "
                                 "maxThreadCount must be positive or equal to AUTO");
    m_maxThreadCount = maxThreadCount;
}

} // namespace Fiber
