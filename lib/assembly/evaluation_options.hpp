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


#ifndef bempp_evaluation_options_hpp
#define bempp_evaluation_options_hpp

#include "../common/common.hpp"

#include "../common/deprecated.hpp"
#include "../fiber/opencl_options.hpp"
#include "../fiber/parallelization_options.hpp"

namespace Bempp
{

using Fiber::OpenClOptions;
using Fiber::ParallelizationOptions;

/** \ingroup potential_operators
 *  \brief Options controlling evaluation of potentials.
 */
class EvaluationOptions
{
public:
    /** \brief Constructor. */
    EvaluationOptions();

    enum { AUTO = -1 };

    // Temporarily removed (OpenCl support is broken).
    // void enableOpenCl(const OpenClOptions& openClOptions);
    // void disableOpenCl();

    /** \brief Set the maximum number of threads used during evaluation of potentials.
     *
     *  \p maxThreadCount must be a positive number or \p AUTO. In the latter
     *  case the number of threads is determined automatically. */
    void setMaxThreadCount(int maxThreadCount);

    /** \brief Set the maximum number of threads used during evaluation of potentials.
     *
     *  \deprecated Use setMaxThreadCount() instead. */
    BEMPP_DEPRECATED void switchToTbb(int maxThreadCount = AUTO);

    /** \brief Return current parallelization options. */
    const ParallelizationOptions& parallelizationOptions() const;

private:
    ParallelizationOptions m_parallelizationOptions;
};

} // namespace Bempp

#endif
