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

#ifndef fiber_cuda_laplace_3d_double_layer_potential_kernel_functor_hpp
#define fiber_cuda_laplace_3d_double_layer_potential_kernel_functor_hpp

#include "cuda_kernel_functor.hpp"

#include "../common/scalar_traits.hpp"

namespace Fiber {

template <typename ValueType_>
class CudaLaplace3dDoubleLayerPotentialKernelFunctor :
    public CudaKernelFunctor<ValueType_> {

public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  __device__ virtual void evaluate(CoordinateType testPointCoo[3],
                                   CoordinateType trialPointCoo[3],
                                   CoordinateType testElemNormal[3],
                                   CoordinateType trialElemNormal[3],
                                   ValueType &result) const {

    CoordinateType numeratorSum = 0., distanceSq = 0.;
#pragma unroll
    for (int coordIndex = 0; coordIndex < 3; ++coordIndex) {
      CoordinateType diff =
          trialPointCoo[coordIndex] - testPointCoo[coordIndex];
      distanceSq += diff * diff;
      numeratorSum += diff * trialElemNormal[coordIndex];
    }
    CoordinateType distance = sqrt(distanceSq);
    result = -numeratorSum / (static_cast<CoordinateType>(4. * M_PI) *
                              distance * distanceSq);
  }
};

} // namespace Fiber

#endif
