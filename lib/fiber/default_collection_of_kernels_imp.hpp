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

#include <boost/utility/enable_if.hpp>
#include <stdexcept>

#define FIBER_HAS_MEM_FUNC(func, name)                                        \
    template<typename T, typename Sign>                                 \
    struct name {                                                       \
        typedef char yes[1];                                            \
        typedef char no [2];                                            \
        template <typename U, U> struct type_check;                     \
        template <typename _1> static yes &chk(type_check<Sign, &_1::func> *); \
        template <typename   > static no  &chk(...);                    \
        static bool const value = sizeof(chk<T>(0)) == sizeof(yes);     \
    }

namespace Fiber
{

FIBER_HAS_MEM_FUNC(estimateRelativeScale, hasEstimateRelativeScale);

//template <class Type>
//class TypeHasEstimateRelativeScale
//{
//    // This type won't compile if the second template parameter isn't of type T,
//    // so I can put a function pointer type in the first parameter and the function
//    // itself in the second thus checking that the function has a specific signature.
//    template <typename T, T> struct TypeCheck;

//    typedef char Yes;
//    typedef long No;

//    // A helper struct to hold the declaration of the function pointer.
//    // Change it if the function signature changes.
//    template <typename T> struct EstimateRelativeScale
//    {
//        typedef typename Type::CoordinateType (T::*fptr)(typename Type::CoordinateType) const;
//    };

//    template <typename T> static Yes HasEstimateRelativeScale(TypeCheck< typename EstimateRelativeScale<T>::fptr, &T::estimateRelativeScale >*);
//    template <typename T> static No  HasEstimateRelativeScale(...);

//public:
//    static bool const value = (sizeof(HasEstimateRelativeScale<Type>(0)) == sizeof(Yes));
//};


template<typename Functor>
typename boost::enable_if<hasEstimateRelativeScale<Functor,
                          typename Functor::CoordinateType(Functor::*)(typename Functor::CoordinateType) const>,
                          typename Functor::CoordinateType>::type
estimateRelativeScaleInternal(const Functor& functor,
                              typename Functor::CoordinateType distance)
{
    return functor.estimateRelativeScale(distance);
}

template<typename Functor>
typename boost::disable_if<hasEstimateRelativeScale<Functor,
                           typename Functor::CoordinateType(Functor::*)(typename Functor::CoordinateType) const>,
                           typename Functor::CoordinateType>::type
estimateRelativeScaleInternal(const Functor& functor,
                              typename Functor::CoordinateType distance)
{
    return 1.;
}

//template<typename Functor>
//typename boost::enable_if<TypeHasEstimateRelativeScale<Functor>,
//                          typename Functor::CoordinateType>::type
//estimateRelativeScaleInternal(const Functor& functor,
//                              typename Functor::CoordinateType distance)
//{
//   // std::cout << "true impl" << std::endl;
//   return functor.estimateRelativeScale(distance);
//}

//template<typename Functor>
//typename boost::disable_if<TypeHasEstimateRelativeScale<Functor>,
//                           typename Functor::CoordinateType>::type
//estimateRelativeScaleInternal(const Functor& functor,
//                              typename Functor::CoordinateType distance)
//{  // std::cout << "fallback" << std::endl;
//   return 1.;
//}

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

template <typename Functor>
void DefaultCollectionOfKernels<Functor>::evaluateOnGrid(
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

#pragma ivdep
    for (size_t trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
        for (size_t testIndex = 0; testIndex < testPointCount; ++testIndex)
            m_functor.evaluate(testGeomData.const_slice(testIndex),
                               trialGeomData.const_slice(trialIndex),
                               result.slice(testIndex, trialIndex).self());
}

template <typename Functor>
std::pair<const char*, int>
DefaultCollectionOfKernels<Functor>::evaluateClCode() const {
    throw std::runtime_error("DefaultCollectionOfKernels::evaluateClCode(): "
                             "not implemented yet");
}


template <typename Functor>
typename DefaultCollectionOfKernels<Functor>::CoordinateType
DefaultCollectionOfKernels<Functor>::
estimateRelativeScale(CoordinateType distance) const
{
    return estimateRelativeScaleInternal(m_functor, distance);
}

} // namespace Fiber

#endif
