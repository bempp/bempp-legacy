#include "evaluation_options.hpp"

#include <stdexcept>

namespace Bempp
{

EvaluationOptions::EvaluationOptions() :
    m_parallelism(TBB), m_maxThreadCount(AUTO)
{
    m_openClOptions.useOpenCl = false;
}

void EvaluationOptions::switchToOpenCl(const OpenClOptions& openClOptions)
{
    m_parallelism = OPEN_CL;
    m_openClOptions = openClOptions;
    m_openClOptions.useOpenCl = true;
}

void EvaluationOptions::switchToTbb(int maxThreadCount)
{
    m_parallelism = TBB;
    m_openClOptions.useOpenCl = false;
    if (maxThreadCount <= 0 && maxThreadCount != AUTO)
        throw std::runtime_error("EvaluationOptions::switchToTbb(): "
                                 "maxThreadCount must be positive or equal to AUTO");
    m_maxThreadCount = maxThreadCount;
}

} // namespace Bempp

