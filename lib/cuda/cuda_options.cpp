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

namespace Bempp {

CudaOptions::CudaOptions()
    : m_precision("double"), m_elementDataCachingEnabled(true),
      m_kernelDataCachingEnabled(true), m_devices({0,1}), m_quadOrder(4),
      m_blockSize(512), m_chunkElemPairCount(AUTO) {
}

void CudaOptions::setPrecision(std::string precision) {
  if (precision == "double" || precision == "single")
    m_precision = precision;
  else
    throw std::runtime_error(
        "CudaOptions::setPrecision(): "
        "precision must be ""double"" or ""single""");
}

std::string CudaOptions::precision() const { return m_precision; }

void CudaOptions::enableElementDataCaching() {
  m_elementDataCachingEnabled = true;
}

void CudaOptions::disableElementDataCaching() {
  m_elementDataCachingEnabled = false;
}

bool CudaOptions::isElementDataCachingEnabled() const {
  return m_elementDataCachingEnabled;
}

void CudaOptions::enableKernelDataCaching() {
  m_kernelDataCachingEnabled = true;
}

void CudaOptions::disableKernelDataCaching() {
  m_kernelDataCachingEnabled = false;
}

bool CudaOptions::isKernelDataCachingEnabled() const {
  return m_kernelDataCachingEnabled;
}

void CudaOptions::setDevices(std::vector<int> deviceIds) {
  m_devices = deviceIds;
}

const std::vector<int>& CudaOptions::devices() const {
  return m_devices;
}

void CudaOptions::setQuadOrder(int quadOrder) {
  if (quadOrder > 0)
    m_quadOrder = quadOrder;
  else
    throw std::runtime_error(
        "CudaOptions::setQuadOrder(): "
        "quadOrder must be positive integer greater zero");
}

int CudaOptions::quadOrder() const { return m_quadOrder; }

void CudaOptions::setBlockSize(int blockSize) {
  if (blockSize > 0)
    m_blockSize = blockSize;
  else
    throw std::runtime_error(
        "CudaOptions::setBlockSize(): "
        "blockSize must be positive integer greater zero");
}

int CudaOptions::blockSize() const { return m_blockSize; }

void CudaOptions::setChunkElemPairCount(int chunkElemPairCount) {
  if (chunkElemPairCount > 0 || chunkElemPairCount == AUTO)
    m_chunkElemPairCount = chunkElemPairCount;
  else
    throw std::runtime_error(
        "CudaOptions::setChunkElemPairCount(): "
        "chunkElemPairCount must be positive integer greater zero or AUTO");
}

int CudaOptions::chunkElemPairCount() const { return m_chunkElemPairCount; }

} // namespace Fiber
