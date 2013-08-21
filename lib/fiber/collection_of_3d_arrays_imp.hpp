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

#include "collection_of_3d_arrays.hpp"

namespace Fiber
{

// CollectionOf3dArrays

template <typename T>
inline CollectionOf3dArrays<T>::CollectionOf3dArrays() :
    m_size(0)
{
}

template <typename T>
inline CollectionOf3dArrays<T>::CollectionOf3dArrays(size_t size) :
    m_size(size), m_arrays(new _3dArray<T>[size])
{
    // should we initialise to 0?
}

template <typename T>
inline void CollectionOf3dArrays<T>::set_size(size_t new_size)
{
    if (new_size == m_size)
        return;
    m_arrays.reset(new_size == 0 ? 0 : new _3dArray<T>[new_size]);
    m_size = new_size;
}

template <typename T>
inline size_t CollectionOf3dArrays<T>::size() const
{
    return m_size;
}

template <typename T>
inline void CollectionOf3dArrays<T>::fill(const T& value)
{
    for (size_t a = 0; a < m_size; ++a)
        std::fill((*this)[a].begin(), (*this)[a].end(), value);
}

template <typename T>
inline _3dArray<T>& CollectionOf3dArrays<T>::operator[](size_t index)
{
    return array(index);
}

template <typename T>
inline const _3dArray<T>& CollectionOf3dArrays<T>::operator[](size_t index) const
{
    return array(index);
}

template <typename T>
inline CollectionOf2dSlicesOf3dArrays<T> CollectionOf3dArrays<T>::slice(size_t index2)
{
    return CollectionOf2dSlicesOf3dArrays<T>(*this, index2);
}

template <typename T>
inline CollectionOf2dSlicesOfConst3dArrays<T> CollectionOf3dArrays<T>::const_slice(
        size_t index2) const
{
    return CollectionOf2dSlicesOfConst3dArrays<T>(*this, index2);
}

template <typename T>
inline CollectionOf1dSlicesOf3dArrays<T> CollectionOf3dArrays<T>::slice(size_t index1, size_t index2)
{
    return CollectionOf1dSlicesOf3dArrays<T>(*this, index1, index2);
}

template <typename T>
inline CollectionOf1dSlicesOfConst3dArrays<T> CollectionOf3dArrays<T>::const_slice(
        size_t index1, size_t index2) const
{
    return CollectionOf1dSlicesOfConst3dArrays<T>(*this, index1, index2);
}

template <typename T>
inline _3dArray<T>& CollectionOf3dArrays<T>::array(size_t i)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_array_index(i);
#endif
    return m_arrays[i];
}

template <typename T>
inline const _3dArray<T>& CollectionOf3dArrays<T>::array(size_t i) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_array_index(i);
#endif
    return m_arrays[i];
}

template <typename T>
inline void CollectionOf3dArrays<T>::check_array_index(size_t array_index) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (m_size <= array_index)
        throw std::invalid_argument("Invalid array index");
#endif
}

// CollectionOf2dSlicesOf3dArrays

template <typename T>
inline CollectionOf2dSlicesOf3dArrays<T>::
CollectionOf2dSlicesOf3dArrays(CollectionOf3dArrays<T>& collection, size_t index2) :
    m_collection(collection), m_index2(index2)
{}

template <typename T>
inline CollectionOf2dSlicesOf3dArrays<T>& CollectionOf2dSlicesOf3dArrays<T>::
self() {
    return *this;
}

template <typename T>
inline _2dSliceOfConst3dArray<T> CollectionOf2dSlicesOf3dArrays<T>::
operator[](size_t index) const {
    return _2dSliceOfConst3dArray<T>(m_collection[index], m_index2);
}

template <typename T>
inline _2dSliceOf3dArray<T> CollectionOf2dSlicesOf3dArrays<T>::
operator[](size_t index) {
    return _2dSliceOf3dArray<T>(m_collection[index], m_index2);
}

template <typename T>
inline size_t CollectionOf2dSlicesOf3dArrays<T>::size() const {
    return m_collection.size();
}

// CollectionOf2dSlicesOfConst3dArrays

template <typename T>
inline CollectionOf2dSlicesOfConst3dArrays<T>::
CollectionOf2dSlicesOfConst3dArrays(const CollectionOf3dArrays<T>& collection,
                                    size_t index2) :
        m_collection(collection), m_index2(index2)
{}

template <typename T>
inline _2dSliceOfConst3dArray<T> CollectionOf2dSlicesOfConst3dArrays<T>::
operator[](size_t index) const {
    return _2dSliceOfConst3dArray<T>(m_collection[index], m_index2);
}

template <typename T>
inline size_t CollectionOf2dSlicesOfConst3dArrays<T>::size() const {
    return m_collection.size();
}

// CollectionOf1dSlicesOf3dArrays

template <typename T>
inline CollectionOf1dSlicesOf3dArrays<T>::
CollectionOf1dSlicesOf3dArrays(CollectionOf3dArrays<T>& collection,
                               size_t index1, size_t index2) :
    m_collection(collection), m_index1(index1), m_index2(index2)
{}

template <typename T>
inline CollectionOf1dSlicesOf3dArrays<T>& CollectionOf1dSlicesOf3dArrays<T>::
self() {
    return *this;
}

template <typename T>
inline _1dSliceOfConst3dArray<T> CollectionOf1dSlicesOf3dArrays<T>::
operator[](size_t index) const {
    return _1dSliceOfConst3dArray<T>(m_collection[index], m_index1, m_index2);
}

template <typename T>
inline _1dSliceOf3dArray<T> CollectionOf1dSlicesOf3dArrays<T>::
operator[](size_t index) {
    return _1dSliceOf3dArray<T>(m_collection[index], m_index1, m_index2);
}

template <typename T>
inline size_t CollectionOf1dSlicesOf3dArrays<T>::size() const {
    return m_collection.size();
}

// CollectionOf1dSlicesOfConst3dArrays

template <typename T>
inline CollectionOf1dSlicesOfConst3dArrays<T>::
CollectionOf1dSlicesOfConst3dArrays(const CollectionOf3dArrays<T>& collection,
                                    size_t index1, size_t index2) :
        m_collection(collection), m_index1(index1), m_index2(index2)
{}

template <typename T>
inline _1dSliceOfConst3dArray<T> CollectionOf1dSlicesOfConst3dArrays<T>::
operator[](size_t index) const {
    return _1dSliceOfConst3dArray<T>(m_collection[index], m_index1, m_index2);
}

template <typename T>
inline size_t CollectionOf1dSlicesOfConst3dArrays<T>::size() const {
    return m_collection.size();
}

} // namespace Fiber
