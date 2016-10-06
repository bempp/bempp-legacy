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

#ifndef fiber_parallelization_options_hpp
#define fiber_parallelization_options_hpp

#include "../common/common.hpp"

#include "opencl_options.hpp"

namespace Fiber {

/** \brief Parallel operation settings. */
class ParallelizationOptions {
public:
  enum { AUTO = -1 };

  /** \brief Constructor. */
  ParallelizationOptions();

  /** \brief Enable OpenCL-based calculations (currently broken). */
  void enableOpenCl(const OpenClOptions &openClOptions);
  /** \brief Disable OpenCL-based calculations. */
  void disableOpenCl();
  /** \brief Return whether OpenCL-based calculations are enabled. */
  bool isOpenClEnabled() const;
  /** \brief Return current settings controlling operation of the GPU. */
  const OpenClOptions &openClOptions() const;

  /** \brief Set the maximum number of threads used during the assembly.
   *
   *  \p maxThreadCount must be a positive number or \p AUTO. In the latter
   *  case the number of threads is determined automatically by Intel
   *  Threading Building Blocks.*/
  void setMaxThreadCount(int maxThreadCount = AUTO);

  /** \brief Return the maximum number of thread used during the assembly.
   *
   *  The returned value can be a positive number or \p AUTO. In the latter
   *  case the number of threads is determined automatically by
   *  Intel Threading Building Blocks. */
  int maxThreadCount() const;

private:
  bool m_openClEnabled;
  OpenClOptions m_openClOptions;
  int m_maxThreadCount;
};

} // namespace Fiber

#endif
