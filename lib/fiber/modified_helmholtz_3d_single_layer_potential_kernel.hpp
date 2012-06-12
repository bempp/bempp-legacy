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

#ifndef fiber_modified_helmholtz_3d_single_layer_potential_kernel_hpp
#define fiber_modified_helmholtz_3d_single_layer_potential_kernel_hpp

#include "../common/common.hpp"

#include "kernel.hpp"

namespace Fiber
{

/** \ingroup modified_helmholtz_3d
 *  \ingroup fiber
 *  \brief Single-layer-potential kernel for the modified Helmholtz equation in 3D.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. Note that setting \p ValueType to a real
 *  type implies that the wave number will also be purely real.
 *
 *  \see modified_helmholtz_3d
 */
template <typename ValueType>
class ModifiedHelmholtz3dSingleLayerPotentialKernel : public Kernel<ValueType>
{
public:
    typedef typename Kernel<ValueType>::CoordinateType CoordinateType;

    /** \brief Construct the kernel, setting the wave number to \p k.
     *
     * See \ref modified_helmholtz_3d for the definition of the wave number. */
    explicit ModifiedHelmholtz3dSingleLayerPotentialKernel(ValueType k);

    virtual int worldDimension() const { return 3; }
    virtual int domainDimension() const { return 1; }
    virtual int codomainDimension() const { return 1; }

    virtual void addGeometricalDependencies(size_t& testGeomDeps,
                                            size_t& trialGeomDeps) const;

    /** \brief Return the current value of wave number.
     *
     * See \ref modified_helmholtz_3d for the definition of the wave number. */
    ValueType waveNumber() const { return m_waveNumber; }

    /** \brief Set the wave number to \p k.
     *
     * See \ref modified_helmholtz_3d for the definition of the wave number. */
    void setWaveNumber(ValueType k) { m_waveNumber = k; }

    virtual void evaluateAtPointPairs(const GeometricalData<CoordinateType>& testGeomData,
                                      const GeometricalData<CoordinateType>& trialGeomData,
                                      arma::Cube<ValueType>& result) const;

    virtual void evaluateOnGrid(const GeometricalData<CoordinateType>& testGeomData,
                                const GeometricalData<CoordinateType>& trialGeomData,
                                Array4d<ValueType>& result) const;

    virtual std::pair<const char*,int> evaluateClCode () const;

private:
    ValueType evaluateAtPointPair(const arma::Col<CoordinateType>& testPoint,
                                  const arma::Col<CoordinateType>& trialPoint) const;

    ValueType m_waveNumber;
};

} // namespace fiber

#endif 
