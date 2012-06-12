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

#ifndef fiber_parallelisation_options_hpp
#define fiber_parallelisation_options_hpp

#include "../common/common.hpp"

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
