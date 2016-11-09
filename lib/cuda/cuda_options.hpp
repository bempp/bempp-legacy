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

#ifndef fiber_cuda_options_hpp
#define fiber_cuda_options_hpp

#include "../common/common.hpp"

#include <vector>

namespace Fiber {

/** \brief CUDA operation settings. */
class CudaOptions {
public:
  enum { AUTO = 1 };

  /** \brief Constructor. */
  CudaOptions();

  /** \brief Enable element data caching on the device. */
  void enableElemDataCaching();
  /** \brief Disable element data caching on the device. */
  void disableElemDataCaching();
  /** \brief Return whether element data caching on the device is enabled. */
  bool isElemDataCachingEnabled() const;

  /** \brief Set the number of concurrent operation streams on the device.
   *
   *  \p streamCount must be a positive number or \p AUTO. In the latter
   *  case only the default stream is used.*/
  void setStreamCount(int streamCount = AUTO);

  /** \brief Return the number of concurrent operation streams on the device.
   *
   *  The returned value can be a positive number or \p AUTO. In the latter
   *  case only the default stream is used.*/
  int streamCount() const;

  /** \brief Set the devices used during the assembly.
   *
   *  \p deviceIds must be a vector of positive numbers including 0. If no
   *  devices are specified, only device 0 is used by default.*/
  void setDevices(std::vector<int> deviceIds = {0});

  /** \brief Return the devices used during the assembly.
   *
   *  The returned vector can hold positive numbers including 0.  If no
   *  devices are specified, only device 0 is used by default.*/
  const std::vector<int>& devices() const;

private:
  bool m_elemDataCachingEnabled;
  int m_streamCount;
  std::vector<int> m_devices;
};

}

#endif
