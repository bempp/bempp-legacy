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

#include "parallelization_options.hpp"

#include <stdexcept>

namespace Fiber {

ParallelizationOptions::ParallelizationOptions()
    : m_openClEnabled(false), m_maxThreadCount(AUTO) {
  m_openClOptions.useOpenCl = false;
}

void ParallelizationOptions::enableOpenCl(const OpenClOptions &openClOptions) {
  m_openClEnabled = true;
  m_openClOptions = openClOptions;
  m_openClOptions.useOpenCl = true;
}

void ParallelizationOptions::disableOpenCl() {
  m_openClEnabled = false;
  m_openClOptions.useOpenCl = false;
}

bool ParallelizationOptions::isOpenClEnabled() const { return m_openClEnabled; }

const OpenClOptions &ParallelizationOptions::openClOptions() const {
  return m_openClOptions;
}

//void ParallelizationOptions::enableCUDA(const CUDAOptions &CUDAOptions) {
//  m_CUDAEnabled = true;
//  m_CUDAOptions = CUDAOptions;
//}
//
//void ParallelizationOptions::disableCUDA() {
//  m_CUDAEnabled = false;
//}
//
//bool ParallelizationOptions::isCUDAEnabled() const { return m_CUDAEnabled; }
//
//const CUDAOptions &ParallelizationOptions::CUDAOptions() const {
//  return m_CUDAOptions;
//}

void ParallelizationOptions::setMaxThreadCount(int maxThreadCount) {
  if (maxThreadCount <= 0 && maxThreadCount != AUTO)
    throw std::runtime_error(
        "ParallelizationOptions::switchToTbb(): "
        "maxThreadCount must be positive or equal to AUTO");
  m_maxThreadCount = maxThreadCount;
}

int ParallelizationOptions::maxThreadCount() const { return m_maxThreadCount; }

} // namespace Fiber
