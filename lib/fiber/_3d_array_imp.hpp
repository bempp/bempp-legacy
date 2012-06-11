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

#include <stdexcept>

namespace Fiber
{

// _3dArray

template <typename T>
inline _3dArray<T>::_3dArray()
{
    m_extents[0] = 0;
    m_extents[1] = 0;
    m_extents[2] = 0;
    m_storage = 0;
    m_owns = false;
}

template <typename T>
inline _3dArray<T>::_3dArray(int extent0, int extent1, int extent2)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_storage = new T[extent0 * extent1 * extent2];
    m_owns = true;
}

template <typename T>
inline _3dArray<T>::_3dArray(int extent0, int extent1, int extent2, T* data)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_storage = data;
    m_owns = false;
}

template <typename T>
inline _3dArray<T>::~_3dArray()
{
    if (m_owns && m_storage)
        delete[] m_storage;
}

template <typename T>
inline T& _3dArray<T>::operator()(int index0, int index1, int index2)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1, index2);
#endif
    return m_storage[
            index0 +
            m_extents[0] * index1 +
            m_extents[0] * m_extents[1] * index2];
}

template <typename T>
inline const T& _3dArray<T>::operator()(int index0, int index1, int index2) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1, index2);
#endif
    return m_storage[
            index0 +
            m_extents[0] * index1 +
            m_extents[0] * m_extents[1] * index2];
}

template <typename T>
inline int _3dArray<T>::extent(int dimension) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_dimension(dimension);
#endif
    return m_extents[dimension];
}

template <typename T>
inline void _3dArray<T>::set_size(int extent0, int extent1, int extent2)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2);
#endif
    if (extent0 * extent1 * extent2 ==
            m_extents[0] * m_extents[1] * m_extents[2]) {
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
    }
    else {
        if (m_owns)
            delete[] m_storage;
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
        m_storage = new T[extent0 * extent1 * extent2];
        m_owns = true;
    }
}

template <typename T>
inline typename _3dArray<T>::iterator _3dArray<T>::begin()
{
    return m_storage;
}

template <typename T>
inline typename _3dArray<T>::const_iterator _3dArray<T>::begin() const
{
    return m_storage;
}

template <typename T>
inline typename _3dArray<T>::iterator _3dArray<T>::end()
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2];
}

template <typename T>
inline typename _3dArray<T>::const_iterator _3dArray<T>::end() const
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2];
}

#ifdef FIBER_CHECK_ARRAY_BOUNDS
template <typename T>
inline void _3dArray<T>::check_dimension(int dimension) const
{
    if (dimension < 0 || 2 < dimension)
        throw std::invalid_argument("Invalid dimension");
}

template <typename T>
inline void _3dArray<T>::check_extents(int extent0, int extent1, int extent2) const
{
    if (extent0 <= 0 || extent1 <= 0 || extent2 <= 0)
        throw std::length_error("Invalid extent");
}

template <typename T>
inline void _3dArray<T>::check_indices(int index0, int index1, int index2) const
{
    if (index0 < 0 || m_extents[0] <= index0 ||
        index1 < 0 || m_extents[1] <= index1 ||
        index2 < 0 || m_extents[2] <= index2)
        throw std::out_of_range("Invalid index");
}
#endif // FIBER_CHECK_ARRAY_BOUNDS

// _2dSliceOf3dArray

template <typename T>
inline _2dSliceOf3dArray<T>::_2dSliceOf3dArray(_3dArray<T>& array, int index2) :
    m_array(array), m_index2(index2)
{}

template <typename T>
inline _2dSliceOf3dArray<T>& _2dSliceOf3dArray<T>::self() {
    return *this;
}

template <typename T>
inline const T& _2dSliceOf3dArray<T>::operator()(int index0, int index1) const {
    return m_array(index0, index1, m_index2);
}

template <typename T>
inline T& _2dSliceOf3dArray<T>::operator()(int index0, int index1) {
    return m_array(index0, index1, m_index2);
}

template <typename T>
inline int _2dSliceOf3dArray<T>::extent(int dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _2dSliceOf3dArray<T>::check_dimension(int dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (dimension < 0 || 1 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _2dSliceOfConst3dArray

template <typename T>
inline _2dSliceOfConst3dArray<T>::_2dSliceOfConst3dArray(
        const _3dArray<T>& array, int index2) :
    m_array(array), m_index2(index2)
{}

template <typename T>
inline const T& _2dSliceOfConst3dArray<T>::operator()(int index0, int index1) const {
    return m_array(index0, index1, m_index2);
}

template <typename T>
inline int _2dSliceOfConst3dArray<T>::extent(int dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _2dSliceOfConst3dArray<T>::check_dimension(int dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (dimension < 0 || 1 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _1dSliceOf3dArray

template <typename T>
inline _1dSliceOf3dArray<T>::_1dSliceOf3dArray(
        _3dArray<T>& array, int index1, int index2) :
    m_array(array), m_index1(index1), m_index2(index2)
{}

template <typename T>
inline _1dSliceOf3dArray<T>& _1dSliceOf3dArray<T>::self() {
    return *this;
}

template <typename T>
inline const T& _1dSliceOf3dArray<T>::operator()(int index0) const {
    return m_array(index0, m_index1, m_index2);
}

template <typename T>
inline T& _1dSliceOf3dArray<T>::operator()(int index0) {
    return m_array(index0, m_index1, m_index2);
}

template <typename T>
inline int _1dSliceOf3dArray<T>::extent(int dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _1dSliceOf3dArray<T>::check_dimension(int dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (dimension < 0 || 0 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _1dSliceOfConst3dArray

template <typename T>
inline _1dSliceOfConst3dArray<T>::_1dSliceOfConst3dArray(
        const _3dArray<T>& array, int index1, int index2) :
    m_array(array), m_index1(index1), m_index2(index2)
{}

template <typename T>
inline const T& _1dSliceOfConst3dArray<T>::operator()(int index0) const {
    return m_array(index0, m_index1, m_index2);
}

template <typename T>
inline int _1dSliceOfConst3dArray<T>::extent(int dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _1dSliceOfConst3dArray<T>::check_dimension(int dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (dimension < 0 || 0 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

} // namespace Fiber
