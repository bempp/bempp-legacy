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

#ifndef fiber_modified_maxwell_3d_far_field_single_layer_potential_operator_kernel_functor_hpp
#define fiber_modified_maxwell_3d_far_field_single_layer_potential_operator_kernel_functor_hpp

#include "../common/common.hpp"

#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

namespace Fiber
{

/** \ingroup modified_maxwell_3d
 *  \ingroup functors 
 *  \brief Kernel collection functor for the far-field single-layer potential
 *  operator of the modified Maxwell equations in 3D.
 *
 *  The functor evaluates two kernels: the Green's function of the modified
 *  Helmholtz equation, multiplied by m_waveNumber, and the gradient of this
 *  Green's function with respect to the test coordinate, divided by
 *  m_waveNumber.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. Note that setting \p ValueType to a real
 *  type implies that the wave number will also be purely real.
 *
 *  \see modified_maxwell_3d
 */
template <typename ValueType_>
class ModifiedMaxwell3dFarFieldSingleLayerPotentialOperatorKernelFunctor
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    explicit ModifiedMaxwell3dFarFieldSingleLayerPotentialOperatorKernelFunctor(
            ValueType waveNumber) :
        m_waveNumber(waveNumber)
    {}

    int kernelCount() const { return 2; }
    int kernelRowCount(int kernelIndex) const { return kernelIndex == 0 ? 1 : 3; }
    int kernelColCount(int kernelIndex) const { return 1; }

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const {
        testGeomDeps |= GLOBALS;
        trialGeomDeps |= GLOBALS;
    }

    ValueType waveNumber() const { return m_waveNumber; }

    template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
    void evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            CollectionOf2dSlicesOfNdArrays<ValueType>& result) const {
        const int coordCount = 3;

        CoordinateType x_y = 0;
        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
            x_y += testGeomData.global(coordIndex) *
                    trialGeomData.global(coordIndex);
        const ValueType commonFactor =
                static_cast<ValueType>(1.0 / (4.0 * M_PI)) *
                exp(m_waveNumber * x_y);
        result[0](0, 0) = m_waveNumber * commonFactor;
        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
            result[1](coordIndex, 0) = -testGeomData.global(coordIndex) *
                    commonFactor;
    }

private:
    /** \cond PRIVATE */
    ValueType m_waveNumber;
    /** \endcond */
};

} // namespace Fiber

#endif
