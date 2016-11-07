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

#ifndef fiber_laplace_3d_single_layer_potential_kernel_functor_hpp
#define fiber_laplace_3d_single_layer_potential_kernel_functor_hpp

#include "../common/common.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../cuda/cuda_laplace_3d_single_layer_potential_kernel_functor.hpp"

#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

namespace Fiber {

/** \ingroup laplace_3d
 *  \ingroup functors
 *  \brief Single-layer-potential kernel functor for the Laplace equation in 3D.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>.
 *
 *  \see laplace_3d
 */

template <typename ValueType_>
class Laplace3dSingleLayerPotentialKernelFunctor {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  shared_ptr<CudaLaplace3dSingleLayerPotentialKernelFunctor<ValueType>>
  cudaFunctor() {
    m_cudaFunctor = boost::make_shared<
        CudaLaplace3dSingleLayerPotentialKernelFunctor<ValueType>>();
    return m_cudaFunctor;
  }

  int kernelCount() const { return 1; }
  int kernelRowCount(int /* kernelIndex */) const { return 1; }
  int kernelColCount(int /* kernelIndex */) const { return 1; }

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS;
  }

  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluate(const ConstGeometricalDataSlice<CoordinateType> &testGeomData,
                const ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    const int coordCount = 3;
    assert(testGeomData.dimWorld() == coordCount);
    assert(result.size() == 1);

    ValueType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      ValueType diff =
          testGeomData.global(coordIndex) - trialGeomData.global(coordIndex);
      sum += diff * diff;
    }
    result[0](0, 0) = static_cast<CoordinateType>(1. / (4. * M_PI)) / sqrt(sum);
  }

private:
  shared_ptr<CudaLaplace3dSingleLayerPotentialKernelFunctor<ValueType>>
  m_cudaFunctor;
};

} // namespace Fiber

#endif
