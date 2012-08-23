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

#ifndef fiber_default_collection_of_kernels_imp_hpp
#define fiber_default_collection_of_kernels_imp_hpp

#include "default_collection_of_kernels.hpp"

#include "collection_of_3d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "geometrical_data.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename Functor>
void DefaultCollectionOfKernels<Functor>::addGeometricalDependencies(
        size_t& testGeomDeps, size_t& trialGeomDeps) const
{
    m_functor.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
}

template <typename Functor>
void DefaultCollectionOfKernels<Functor>::evaluateAtPointPairs(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<ValueType>& result) const
{
    assert(testGeomData.pointCount() == trialGeomData.pointCount());

    const size_t pointCount = testGeomData.pointCount();
    const size_t kernelCount = m_functor.kernelCount();
    result.set_size(kernelCount);
    for (size_t k = 0; k < kernelCount; ++k)
        result[k].set_size(m_functor.kernelRowCount(k),
                           m_functor.kernelColCount(k),
                           pointCount);

    for (size_t p = 0; p < pointCount; ++p)
        m_functor.evaluate(testGeomData.const_slice(p),
                           trialGeomData.const_slice(p),
                           result.slice(p).self());
}

template <typename ValueType>
void DefaultCollectionOfKernels<ValueType>::evaluateOnGrid(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf4dArrays<ValueType>& result) const
{
    const size_t testPointCount = testGeomData.pointCount();
    const size_t trialPointCount = trialGeomData.pointCount();
    const size_t kernelCount = m_functor.kernelCount();
    result.set_size(kernelCount);
    for (size_t k = 0; k < kernelCount; ++k)
        result[k].set_size(m_functor.kernelRowCount(k),
                           m_functor.kernelColCount(k),
                           testPointCount,
                           trialPointCount);

    for (size_t trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (size_t testIndex = 0; testIndex < testPointCount; ++testIndex)
            m_functor.evaluate(testGeomData.const_slice(testIndex),
                               trialGeomData.const_slice(trialIndex),
                               result.slice(testIndex, trialIndex).self());
}

template <typename ValueType>
std::pair<const char*, int>
DefaultCollectionOfKernels<ValueType>::evaluateClCode() const {
    throw std::runtime_error("DefaultCollectionOfKernels::evaluateClCode(): "
                             "not implemented yet");
}

} // namespace Fiber

#endif
