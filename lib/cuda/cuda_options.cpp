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

#include "cuda_options.hpp"

#include <stdexcept>

namespace Fiber {

CudaOptions::CudaOptions()
    : m_elemDataCachingEnabled(false), m_streamCount(AUTO) {
  m_devices.resize(1);
  m_devices[0] = 0;
}

void CudaOptions::enableElemDataCaching() {
  m_elemDataCachingEnabled = true;
}

void CudaOptions::disableElemDataCaching() {
  m_elemDataCachingEnabled = false;
}

bool CudaOptions::isElemDataCachingEnabled() const {
  return m_elemDataCachingEnabled;
}

void CudaOptions::setStreamCount(int streamCount) {
  if (streamCount < 1)
    throw std::runtime_error(
        "CudaOptions::setStreamCount(): "
        "streamCount must be positive integer greater zero");
  m_streamCount = streamCount;
}

int CudaOptions::streamCount() const { return m_streamCount; }

void CudaOptions::setDevices(std::vector<int> deviceIds) {
  m_devices = deviceIds;
}

const std::vector<int>& CudaOptions::devices() const {
  return m_devices;
}

} // namespace Fiber
