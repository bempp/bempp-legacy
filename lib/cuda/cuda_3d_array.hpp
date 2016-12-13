// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_cuda_3d_array_hpp
#define fiber_cuda_3d_array_hpp

#include <thrust/device_vector.h>

namespace Fiber {

template <typename T>
class Cuda3dArray {

public:

  Cuda3dArray() : m_extent0(0), m_extent1(0), m_extent2(0) { }

  Cuda3dArray(const size_t extent0, const size_t extent1, const size_t extent2)
              : m_extent0(extent0), m_extent1(extent1), m_extent2(extent2) {

    resize(m_extent0, m_extent1, m_extent2);
  }

  ~Cuda3dArray() { }

  void resize(const size_t newExtent0,
              const size_t newExtent1,
              const size_t newExtent2) {

    m_extent0 = newExtent0;
    m_extent1 = newExtent1;
    m_extent2 = newExtent2;

    m_data.resize(m_extent0 * m_extent1 * m_extent2);
  }

  T operator()(const size_t index0, const size_t index1, const size_t index2) const {

    return m_data[getIndex(index0, index1, index2)];
  }

  T& operator()(const size_t index0, const size_t index1, const size_t index2) {

    return m_data[getIndex(index0, index1, index2)];
  }

private:

  size_t getIndex(const size_t index0, const size_t index1, const size_t index2) const {

    // Index dim2 varies fastest in data array
    return index0 * m_extent1 * m_extent2
         + index1 * m_extent2
         + index2;
  }

  size_t m_extent0, m_extent1, m_extent2;

  thrust::device_vector<T> m_data;
};

} // namespace Fiber

#endif
