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

namespace Bempp {

/** \brief CUDA operation settings. */
class CudaOptions {
public:

  /** \brief Constructor. */
  CudaOptions();

  /** \brief Set the data precision used on the device. */
  void setPrecision(std::string precision);
  /** \brief Return the data precision used on the device. */
  std::string precision() const;

  /** \brief Enable element data caching on the device. */
  void enableElementDataCaching();
  /** \brief Disable element data caching on the device. */
  void disableElementDataCaching();
  /** \brief Return whether element data caching on the device is enabled. */
  bool isElementDataCachingEnabled() const;

  /** \brief Set the number of concurrent operation streams on the device. */
  void setStreamCount(int streamCount);
  /** \brief Return the number of concurrent operation streams on the device. */
  int streamCount() const;

  /** \brief Set the devices used during the assembly. */
  void setDevices(std::vector<int> deviceIds = {0});
  /** \brief Return the devices used during the assembly. */
  const std::vector<int>& devices() const;

  /** \brief Set the order used for numerical quadrature on the device. */
  void setQuadOrder(int quadOrder);
  /** \brief Return the order used for numerical quadrature on the device. */
  int quadOrder() const;

  /** \brief Set the block size used for calculations on the device. */
  void setBlockSize(int blockSize);
  /** \brief Return the block size used for calculations on the device. */
  int blockSize() const;

private:
  std::string m_precision;
  bool m_elementDataCachingEnabled;
  int m_streamCount;
  std::vector<int> m_devices;
  int m_quadOrder;
  int m_blockSize;
};

}

#endif
