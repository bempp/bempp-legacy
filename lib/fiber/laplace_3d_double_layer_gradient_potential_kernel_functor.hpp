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

#ifndef fiber_laplace_3d_double_layer_gradient_potential_kernel_functor_hpp
#define fiber_laplace_3d_double_layer_gradient_potential_kernel_functor_hpp

#include "../common/common.hpp"

#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

namespace Fiber {

/** \ingroup laplace_3d
 *  \ingroup functors
 *  \brief Double-layer-gradient-potential kernel functor for the Laplace
 * equation in 3D.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>.
 *
 *  \see laplace_3d
 */

template <typename ValueType_>
class Laplace3dDoubleLayerGradientPotentialKernelFunctor {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  int kernelCount() const { return 1; }
  int kernelRowCount(int /* kernelIndex */) const { return 3; }
  int kernelColCount(int /* kernelIndex */) const { return 1; }

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS | NORMALS;
  }

  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluate(const ConstGeometricalDataSlice<CoordinateType> &testGeomData,
                const ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    const int coordCount = 3;
    assert(testGeomData.dimWorld() == coordCount);
    assert(result.size() == 3);

    CoordinateType sum = 0;
    CoordinateType diff[coordCount];
    CoordinateType testGlobal[coordCount];
    CoordinateType normal[coordCount];

    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      normal[coordIndex] = trialGeomData.normal(coordIndex);
      testGlobal[coordIndex] = testGeomData.global(coordIndex);
      diff[coordIndex] =
          testGlobal[coordIndex] - trialGeomData.global(coordIndex);
      sum += diff[coordIndex] * diff[coordIndex];
    }
    CoordinateType sqrt_sum = sqrt(sum);

    CoordinateType term1[coordCount];
    CoordinateType term2[coordCount];
    CoordinateType factor = 0;

    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
      factor += normal[coordIndex] * diff[coordIndex] * 3. / sum;

    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      term1[coordIndex] = normal[coordIndex];
      term2[coordIndex] = factor * diff[coordIndex];
    }

    result[0](0, 0) = static_cast<CoordinateType>(1. / (4. * M_PI)) *
                      (term1[0] - term2[0]) / (sum * sqrt_sum);
    result[0](1, 0) = static_cast<CoordinateType>(1. / (4. * M_PI)) *
                      (term1[1] - term2[1]) / (sum * sqrt_sum);
    result[0](2, 0) = static_cast<CoordinateType>(1. / (4. * M_PI)) *
                      (term1[2] - term2[2]) / (sum * sqrt_sum);
  }
};

} // namespace Fiber

#endif
