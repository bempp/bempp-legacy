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

#ifndef fiber_kernel_hpp
#define fiber_kernel_hpp

#include "array_4d.hpp"

#include <armadillo>

namespace Fiber
{

template <typename ValueType> class BasisData;
template <typename ValueType> class GeometricalData;

template <typename ValueType>
class Kernel
{
public:
    virtual int worldDimension() const = 0;
    virtual int domainDimension() const = 0;
    virtual int codomainDimension() const = 0;

    virtual void addGeometricalDependencies(int& testGeomDeps,
                                            int& trialGeomDeps) const = 0;

    virtual void evaluateAtPointPairs(const GeometricalData<ValueType>& trialGeomData,
                                      const GeometricalData<ValueType>& testGeomData,
                                      arma::Cube<ValueType>& result) const = 0;

    virtual void evaluateOnGrid(const GeometricalData<ValueType>& trialGeomData,
                                const GeometricalData<ValueType>& testGeomData,
                                Array4D<ValueType>& result) const = 0;
};

} // namespace Fiber

#endif
