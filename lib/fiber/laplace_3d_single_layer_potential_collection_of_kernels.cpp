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

#include "laplace_3d_single_layer_potential_collection_of_kernels.hpp"

#include "../common/complex_aux.hpp"
#include "collection_of_3d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "explicit_instantiation.hpp"
#include "geometrical_data.hpp"

#include "../common/armadillo_fwd.hpp"
#include <boost/utility/enable_if.hpp>
#include <stdexcept>

namespace Fiber
{

template <typename ValueType>
void Laplace3dSingleLayerPotentialCollectionOfKernels<ValueType>::addGeometricalDependencies(
        size_t& testGeomDeps, size_t& trialGeomDeps) const
{
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS;
}

template <typename ValueType>
void Laplace3dSingleLayerPotentialCollectionOfKernels<ValueType>::evaluateAtPointPairs(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<ValueType>& result) const
{
    assert(testGeomData.pointCount() == trialGeomData.pointCount());

    const size_t pointCount = testGeomData.pointCount();
    const size_t kernelCount = 1;
    result.set_size(kernelCount);
    for (size_t k = 0; k < kernelCount; ++k)
        result[k].set_size(1, // row count
                           1, // column count
                           pointCount);

    ValueType* __restrict result_ = result[0].begin();
    const CoordinateType* __restrict testGlobals = testGeomData.globals.memptr();
    const CoordinateType* __restrict trialGlobals = trialGeomData.globals.memptr();
    // std::cout << (size_t)result_ % 16 << std::endl;
    // std::cout << (size_t)testGlobals % 16 << std::endl;
    // std::cout << (size_t)trialGlobals % 16 << std::endl;
    // __assume_aligned(result_, 16);
    // __assume_aligned(testGlobals, 16);
    // __assume_aligned(trialGlobals, 16);

#pragma ivdep
// #pragma vector aligned
    for (size_t p = 0; p < pointCount; ++p) {
        const int coordCount = 3;
        
        CoordinateType sum = 0;
        for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
        {
            CoordinateType diff =
                testGlobals[coordIndex + coordCount * p] -
                trialGlobals[coordIndex + coordCount * p];
            sum += diff * diff;
        }
        CoordinateType distance = sqrt(sum);
        result_[p] =
            static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance;
    }
}

template <typename ValueType>
void Laplace3dSingleLayerPotentialCollectionOfKernels<ValueType>::evaluateOnGrid(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf4dArrays<ValueType>& result) const
{
    const size_t testPointCount = testGeomData.pointCount();
    const size_t trialPointCount = trialGeomData.pointCount();
    const size_t kernelCount = 1;
    result.set_size(kernelCount);
    for (size_t k = 0; k < kernelCount; ++k)
        result[k].set_size(1, // row count
                           1, // column count
                           testPointCount,
                           trialPointCount);

    ValueType* __restrict result_ = result[0].begin();
    const CoordinateType* __restrict testGlobals = testGeomData.globals.memptr();
    const CoordinateType* __restrict trialGlobals = trialGeomData.globals.memptr();
    // std::cout << (size_t)result_ % 16 << std::endl;
    // std::cout << (size_t)testGlobals % 16 << std::endl;
    // std::cout << (size_t)trialGlobals % 16 << std::endl;
    // __assume_aligned(result_, 16);
    // __assume_aligned(testGlobals, 16);
    // __assume_aligned(trialGlobals, 16);

    const int coordCount = 3;
    // const int four = 4;
    // CoordinateType arma_aligned testGlobals[256];
    // CoordinateType arma_aligned trialGlobals[256];
    // for (size_t p = 0; p < testPointCount; ++p)
    //     for (int c = 0; c < coordCount; ++c)
    //         testGlobals[c + p * four] = testGeomData.globals(c, p);
    // for (size_t p = 0; p < trialPointCount; ++p)
    //     for (int c = 0; c < coordCount; ++c)
    //         trialGlobals[c + p * four] = trialGeomData.globals(c, p);

#pragma ivdep
//#pragma vector aligned
    for (size_t trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
       for (size_t testIndex = 0; testIndex < testPointCount; ++testIndex) {

           CoordinateType sum = 0;
           for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex)
           {
               CoordinateType diff =
                   testGlobals[coordIndex + coordCount * testIndex] -
                   trialGlobals[coordIndex + coordCount * trialIndex];
                   // testGlobals[coordIndex + four * testIndex] -
                   // trialGlobals[coordIndex + four * trialIndex];
               sum += diff * diff;
           }
           CoordinateType distance = sqrt(sum);
           result_[testIndex + testPointCount * trialIndex] =
               static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance;
       }
}

template <typename ValueType>
std::pair<const char*, int>
Laplace3dSingleLayerPotentialCollectionOfKernels<ValueType>::evaluateClCode() const {
    throw std::runtime_error("DefaultCollectionOfKernels::evaluateClCode(): "
                             "not implemented yet");
}


template <typename ValueType>
typename Laplace3dSingleLayerPotentialCollectionOfKernels<ValueType>::CoordinateType
Laplace3dSingleLayerPotentialCollectionOfKernels<ValueType>::
estimateRelativeScale(CoordinateType distance) const
{
    return 1.;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(Laplace3dSingleLayerPotentialCollectionOfKernels);

} // namespace Fiber
