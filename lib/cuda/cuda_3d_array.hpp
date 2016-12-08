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

#include <thrust/host_vector.h>

namespace Fiber {

template <typename T>
class Cuda3dArray {

public:

  Cuda3dArray() : dim1size(0), dim2size(0), dim3size(0) { }

  Cuda3dArray(const unsigned int _dim1size,
              const unsigned int _dim2size,
              const unsigned int _dim3size)
              : dim1size(_dim1size),
                dim2size(_dim2size),
                dim3size(_dim3size) { resize(dim1size, dim2size, dim3size); }

  ~Cuda3dArray() { }

  void resize(const unsigned int _dim1size,
              const unsigned int _dim2size,
              const unsigned int _dim3size) {

    dim1size = _dim1size;
    dim2size = _dim2size;
    dim2size = _dim2size;

    data.resize(dim1size * dim2size * dim3size);
  }

  void clear() { data.clear(); }

  T operator()(const unsigned int _dim1,
               const unsigned int _dim2,
               const unsigned int _dim3) const {

    return data[getIndex(_dim1, _dim2, _dim3)];
  }

  T& operator()(const unsigned int _dim1,
                const unsigned int _dim2,
                const unsigned int _dim3) {

    return data[getIndex(_dim1, _dim2, _dim3)];
  }

private:

  unsigned int dim1size;
  unsigned int dim2size;
  unsigned int dim3size;

  thrust::host_vector<T> data;

  unsigned int getIndex(const unsigned int _dim1,
                        const unsigned int _dim2,
                        const unsigned int _dim3) const {

    // Index dim3 varies fastest in data array
    return _dim1 * dim2size * dim3size
         + _dim2 * dim3size
         + _dim3;
  }
};

} // namespace Fiber

#endif
