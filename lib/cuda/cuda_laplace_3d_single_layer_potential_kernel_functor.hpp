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

#ifndef fiber_cuda_laplace_3d_single_layer_potential_kernel_functor_hpp
#define fiber_cuda_laplace_3d_single_layer_potential_kernel_functor_hpp

#include "../common/scalar_traits.hpp"

namespace Fiber {

template <typename ValueType_>
class CudaLaplace3dSingleLayerPotentialKernelFunctor {

public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

#ifdef __CUDACC__
  __device__ static void evaluate(CoordinateType testGeomData[3],
                                  CoordinateType trialGeomData[3],
                                  ValueType &result) {
    ValueType sum = 0;
    for (int coordIndex = 0; coordIndex < 3; ++coordIndex) {
      ValueType diff = testGeomData[coordIndex] - trialGeomData[coordIndex];
      sum += diff * diff;
    }
    result = static_cast<CoordinateType>(1. / (4. * M_PI)) / sqrt(sum);
  }
#endif
};

} // namespace Fiber

#endif
