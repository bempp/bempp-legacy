#ifndef bempp_evaluation_options_hpp
#define bempp_evaluation_options_hpp

#include "../fiber/opencl_options.hpp"

namespace Bempp
{

using Fiber::OpenClOptions;

class EvaluationOptions
{
public:
    EvaluationOptions();

    enum Mode {
        AUTO = -1,
        NO = 0,
        YES = 1
    };

    /** @}
      @name Parallelism
      @{ */

    enum Parallelism {
        TBB, OPEN_CL
    };

    void switchToOpenCl(const OpenClOptions& openClOptions);
    void switchToTbb(int maxThreadCount = AUTO);

    Parallelism parallelism() const {
        return m_parallelism;
    }

    const OpenClOptions& openClOptions() const {
        return m_openClOptions;
    }

    int maxThreadCount() const {
        return m_maxThreadCount;
    }

private:
    Parallelism m_parallelism;
    OpenClOptions m_openClOptions;
    int m_maxThreadCount;
};

} // namespace Bempp

#endif
