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

#include "collection_of_2d_arrays.hpp"

namespace Fiber
{

// CollectionOf2dArrays

template <typename T>
inline CollectionOf2dArrays<T>::CollectionOf2dArrays() :
    m_size(0)
{
}

template <typename T>
inline CollectionOf2dArrays<T>::CollectionOf2dArrays(int size) :
    m_size(size), m_arrays(new _2dArray<T>[size])
{
    // should we initialise to 0?
}


template <typename T>
inline CollectionOf2dArrays<T>::CollectionOf2dArrays(
        const CollectionOf2dArrays& other) :
    m_size(other.m_size), m_arrays(new _2dArray<T>[other.m_size])
{
    for (int a = 0; a < m_size; ++a)
        m_arrays[a] = other.array(a);
}

template <typename T>
inline CollectionOf2dArrays<T>& CollectionOf2dArrays<T>::operator=(
        const CollectionOf2dArrays& rhs)
{
    if (&rhs != this) {
        set_size(rhs.m_size);
        for (int a = 0; a < m_size; ++a)
            m_arrays[a] = rhs.array(a);
    }
    return *this;
}

template <typename T>
inline void CollectionOf2dArrays<T>::set_size(int new_size)
{
    if (new_size == m_size)
        return;
    m_arrays.reset(new_size == 0 ? 0 : new _2dArray<T>[new_size]);
    m_size = new_size;
}

template <typename T>
inline int CollectionOf2dArrays<T>::size() const
{
    return m_size;
}

template <typename T>
inline void CollectionOf2dArrays<T>::fill(const T& value)
{
    for (int a = 0; a < m_size; ++a)
        std::fill((*this)[a].begin(), (*this)[a].end(), value);
}

template <typename T>
inline _2dArray<T>& CollectionOf2dArrays<T>::operator[](int index)
{
    return array(index);
}

template <typename T>
inline const _2dArray<T>& CollectionOf2dArrays<T>::operator[](int index) const
{
    return array(index);
}

template <typename T>
inline CollectionOf1dSlicesOf2dArrays<T> CollectionOf2dArrays<T>::slice(int index1)
{
    return CollectionOf1dSlicesOf2dArrays<T>(*this, index1);
}

template <typename T>
inline CollectionOf1dSlicesOfConst2dArrays<T> CollectionOf2dArrays<T>::const_slice(
        int index1) const
{
    return CollectionOf1dSlicesOfConst2dArrays<T>(*this, index1);
}

template <typename T>
inline _2dArray<T>& CollectionOf2dArrays<T>::array(int i)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_array_index(i);
#endif
    return m_arrays[i];
}

template <typename T>
inline const _2dArray<T>& CollectionOf2dArrays<T>::array(int i) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_array_index(i);
#endif
    return m_arrays[i];
}

template <typename T>
inline void CollectionOf2dArrays<T>::check_array_index(int array_index) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (array_index < 0 || m_size <= array_index)
        throw std::invalid_argument("Invalid array index");
#endif
}

// CollectionOf1dSlicesOf2dArrays

template <typename T>
inline CollectionOf1dSlicesOf2dArrays<T>::
CollectionOf1dSlicesOf2dArrays(CollectionOf2dArrays<T>& collection,
                               int index1) :
    m_collection(collection), m_index1(index1)
{}

template <typename T>
inline CollectionOf1dSlicesOf2dArrays<T>& CollectionOf1dSlicesOf2dArrays<T>::
self() {
    return *this;
}

template <typename T>
inline _1dSliceOfConst2dArray<T> CollectionOf1dSlicesOf2dArrays<T>::
operator[](int index) const {
    return _1dSliceOfConst2dArray<T>(m_collection[index], m_index1);
}

template <typename T>
inline _1dSliceOf2dArray<T> CollectionOf1dSlicesOf2dArrays<T>::
operator[](int index) {
    return _1dSliceOf2dArray<T>(m_collection[index], m_index1);
}

template <typename T>
inline int CollectionOf1dSlicesOf2dArrays<T>::size() const {
    return m_collection.size();
}

// CollectionOf1dSlicesOfConst2dArrays

template <typename T>
inline CollectionOf1dSlicesOfConst2dArrays<T>::
CollectionOf1dSlicesOfConst2dArrays(const CollectionOf2dArrays<T>& collection,
                                    int index1) :
        m_collection(collection), m_index1(index1)
{}

template <typename T>
inline _1dSliceOfConst2dArray<T> CollectionOf1dSlicesOfConst2dArrays<T>::
operator[](int index) const {
    return _1dSliceOfConst2dArray<T>(m_collection[index], m_index1);
}

template <typename T>
inline int CollectionOf1dSlicesOfConst2dArrays<T>::size() const {
    return m_collection.size();
}

} // namespace Fiber
