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

#include "kernel.hpp"

namespace Fiber
{

template <typename ValueType>
class ModifiedHelmholtz3dSingleLayerPotentialKernel : public Kernel<ValueType>
{
public:
    typedef typename Kernel<ValueType>::CoordinateType CoordinateType;

    virtual int worldDimension() const { return 3; }
    virtual int domainDimension() const { return 1; }
    virtual int codomainDimension() const { return 1; }

    virtual void addGeometricalDependencies(int& testGeomDeps,
                                            int& trialGeomDeps) const;

    virtual ValueType waveNumber() const { return m_waveNumber; }
    virtual void setWaveNumber(ValueType k) { m_waveNumber = k; }

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
