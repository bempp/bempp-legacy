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

#ifndef fiber_collection_of_kernels_hpp
#define fiber_collection_of_kernels_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"

#include <utility>

namespace Fiber
{

template <typename T> class CollectionOf3dArrays;
template <typename T> class CollectionOf4dArrays;
template <typename CoordinateType> class GeometricalData;

template <typename ValueType_>
class CollectionOfKernels
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    virtual ~CollectionOfKernels()
    {}

    virtual void addGeometricalDependencies(
            size_t& testGeomDeps, size_t& trialGeomDeps) const = 0;

    virtual void evaluateAtPointPairs(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            CollectionOf3dArrays<ValueType>& result) const = 0;

    virtual void evaluateOnGrid(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            CollectionOf4dArrays<ValueType>& result) const = 0;

    /**
     * \brief Returns an OpenCL code snippet for kernel evaluation as a string.
     * \note The code snippet must provide device functions devKernevalGrid and
     *    devKernevalPair (for an implementation example, see
     *    CL/laplace_3d_single_layer_potential_kernel.cl)
     * \note The required data must have been pushed to device memory before
     *    invocation.
     */
    virtual std::pair<const char*, int> evaluateClCode() const = 0;
};

} // namespace Fiber

#endif
