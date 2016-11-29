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

#ifndef fiber_modified_helmholtz_3d_hypersingular_kernel_interpolated_functor_hpp
#define fiber_modified_helmholtz_3d_hypersingular_kernel_interpolated_functor_hpp

#include "../common/common.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include "../cuda/cuda_modified_helmholtz_3d_hypersingular_kernel_interpolated_functor.hpp"

#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

#include "modified_helmholtz_3d_single_layer_potential_kernel_interpolated_functor.hpp"

namespace Fiber {

/** \ingroup modified_helmholtz_3d
 *  \ingroup functors
 *  \brief Hypersingular kernel collection functor for the modified Helmholtz
 *  equation in 3D.
 *
 *  The functor evaluates two kernels: the single-layer potential kernel
 *  the single-layer potential kernel multiplied by m_waveNumber**2.
 *
 *  Uses interpolation to speed up calculations.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. Note that setting \p ValueType to a real
 *  type implies that the wave number will also be purely real.
 *
 *  \see modified_helmholtz_3d
 */

template <typename ValueType_>
class ModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  ModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor(
      ValueType waveNumber, CoordinateType maxDist, int interpPtsPerWavelength)
      : m_slpKernel(waveNumber, maxDist, interpPtsPerWavelength) {}

  shared_ptr<const CudaModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor<ValueType>>
  cudaFunctor() {
    m_cudaFunctor = boost::make_shared<
        CudaModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor<ValueType>>();
    return m_cudaFunctor;
  }

  int kernelCount() const { return 2; }
  int kernelRowCount(int /* kernelIndex */) const { return 1; }
  int kernelColCount(int /* kernelIndex */) const { return 1; }

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    m_slpKernel.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  }

  ValueType waveNumber() const { return m_slpKernel.waveNumber(); }

  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluate(const ConstGeometricalDataSlice<CoordinateType> &testGeomData,
                const ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    // This will put the value of the SLP kernel in result[0](0, 0)
    m_slpKernel.evaluate(testGeomData, trialGeomData, result);
    result[1](0, 0) =
        result[0](0, 0) * m_slpKernel.waveNumber() * m_slpKernel.waveNumber();
  }

  CoordinateType estimateRelativeScale(CoordinateType distance) const {
    return m_slpKernel.estimateRelativeScale(distance);
  }

private:
  /** \cond PRIVATE */
  ModifiedHelmholtz3dSingleLayerPotentialKernelInterpolatedFunctor<ValueType>
      m_slpKernel;
  shared_ptr<CudaModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor<ValueType>>
  m_cudaFunctor;
  /** \endcond */
};

} // namespace Fiber

#endif
