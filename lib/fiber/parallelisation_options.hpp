#ifndef fiber_parallelisation_options_hpp
#define fiber_parallelisation_options_hpp

#include "opencl_options.hpp"

namespace Fiber
{

class ParallelisationOptions
{
public:
    enum { AUTO = -1 };
    enum Mode { TBB, OPEN_CL };

    ParallelisationOptions();

    void switchToOpenCl(const OpenClOptions& openClOptions);
    void switchToTbb(int maxThreadCount = AUTO);

    Mode mode() const {
        return m_mode;
    }

    const OpenClOptions& openClOptions() const {
        return m_openClOptions;
    }

    int maxThreadCount() const {
        return m_maxThreadCount;
    }

private:
    Mode m_mode;
    OpenClOptions m_openClOptions;
    int m_maxThreadCount;
};

} // namespace Fiber

#endif
