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

#ifndef fiber_double_layer_potential_3d_kernel_hpp
#define fiber_double_layer_potential_3d_kernel_hpp

#include "kernel.hpp"
#include <armadillo>

namespace Fiber
{

template <typename ValueType>
class DoubleLayerPotential3DKernel : public Kernel<ValueType>
{
public:
    virtual int worldDimension() const { return 3; }
    virtual int domainDimension() const { return 1; }
    virtual int codomainDimension() const { return 1; }

    virtual void addGeometricalDependencies(int& testGeomDeps,
                                            int& trialGeomDeps) const;

    virtual void evaluateAtPointPairs(const GeometricalData<ValueType>& testGeomData,
                                      const GeometricalData<ValueType>& trialGeomData,
                                      arma::Cube<ValueType>& result) const;

    virtual void evaluateOnGrid(const GeometricalData<ValueType>& testGeomData,
                                const GeometricalData<ValueType>& trialGeomData,
                                Array4D<ValueType>& result) const;

    /**
     * \brief Returns an OpenCL code snippet for kernel evaluation as a string.
     * \note The code snippet provides device functions devKernevalGrid and devKernevalPair
     *   (see CL/single_layer_potential_3D_kernel.cl)
     * \note This method is independent of data, unlike the CPU versions, because the
     *   data for multiple elements are pushed to the device separately.
     */
    virtual std::string evaluateClCode () const;

private:
    ValueType evaluateAtPointPair(const arma::Col<ValueType>& testPoint,
                                  const arma::Col<ValueType>& trialPoint,
                                  const arma::Col<ValueType>& trialNormal) const;
};

} // namespace Fiber

#endif
