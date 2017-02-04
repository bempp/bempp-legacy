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

#ifndef fiber_cuda_hpp
#define fiber_cuda_hpp

#include "../fiber/scalar_traits.hpp"

#include <thrust/device_ptr.h>

namespace Fiber {

template <typename ValueType>
struct QuadData {
  unsigned int pointCount;
  thrust::device_ptr<ValueType> weights;
};

template <typename ValueType>
struct BasisFunData {
  unsigned int dofCount;
  thrust::device_ptr<ValueType> values;
  thrust::device_ptr<ValueType> derivatives;
};

template <typename ValueType>
struct GeomShapeFunData {
  thrust::device_ptr<const ValueType> fun0;
  thrust::device_ptr<const ValueType> fun1;
  thrust::device_ptr<const ValueType> fun2;
};

template <typename ValueType>
struct RawGeometryData {
  unsigned int elemCount, vtxCount;
  thrust::device_ptr<const ValueType> vertices;
  thrust::device_ptr<const int> elementCorners;
};

template <typename ValueType>
struct ElemData {
  unsigned int activeElemCount;
  thrust::device_ptr<typename ScalarTraits<ValueType>::RealType> geomData;
  thrust::device_ptr<const typename ScalarTraits<ValueType>::RealType> normals;
  thrust::device_ptr<const typename ScalarTraits<ValueType>::RealType> integrationElements;
  thrust::device_ptr<const ValueType> surfaceCurls;
};

#define cu_verify(x) do {                                                \
    cudaError_t result = x;                                              \
    if (result != cudaSuccess) {                                         \
      fprintf(stderr,"%s:%i: error: cuda function call failed:\n"        \
              "  %s;\nmessage: %s\n",                                    \
              __FILE__,__LINE__,#x,cudaGetErrorString(result));          \
      exit(1);                                                           \
    }                                                                    \
  } while(0)
#define cu_verify_void(x) do {                                           \
    x;                                                                   \
    cudaError_t result = cudaGetLastError();                             \
    if (result != cudaSuccess) {                                         \
      fprintf(stderr,"%s:%i: error: cuda function call failed:\n"        \
              "  %s;\nmessage: %s\n",                                    \
              __FILE__,__LINE__,#x,cudaGetErrorString(result));          \
      exit(1);                                                           \
    }                                                                    \
  } while(0)

} // namespace Fiber

#endif /* fiber_cuda_cuh */
