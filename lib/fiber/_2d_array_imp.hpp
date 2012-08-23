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

template <typename T>
inline _2dArray<T>::_2dArray()
{
    m_storage = 0;
    m_owns = false;
    m_extents[0] = 0;
    m_extents[1] = 0;
}

template <typename T>
inline _2dArray<T>::_2dArray(size_t extent0, size_t extent1)
{
    init_memory(extent0, extent1);
}

template <typename T>
inline _2dArray<T>::_2dArray(size_t extent0, size_t extent1, T* data)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1);
#endif
    m_storage = data;
    m_owns = false;
    m_extents[0] = extent0;
    m_extents[1] = extent1;
}

template <typename T>
inline _2dArray<T>::_2dArray(const _2dArray& other)
{
    init_memory(other.m_extents[0], other.m_extents[1]);
    std::copy(other.begin(), other.end(), m_storage);
}

template <typename T>
inline _2dArray<T>& _2dArray<T>::operator=(const _2dArray& rhs)
{
    if (&rhs != this) {
        set_size(rhs.m_extents[0], rhs.m_extents[1]);
        std::copy(rhs.begin(), rhs.end(), m_storage);
    }
    return *this;
}

template <typename T>
inline _2dArray<T>::~_2dArray()
{
    free_memory();
}

template <typename T>
inline void _2dArray<T>::init_memory(size_t extent0, size_t extent1)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1);
#endif
    m_storage = new T[extent0 * extent1];
    m_owns = true;
    m_extents[0] = extent0;
    m_extents[1] = extent1;
}

template <typename T>
inline void _2dArray<T>::free_memory()
{
    if (m_owns && m_storage)
        delete[] m_storage;
    m_owns = false;
    m_storage = 0;
}

template <typename T>
inline T& _2dArray<T>::operator()(size_t index0, size_t index1)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1);
#endif
    return m_storage[index0 + m_extents[0] * index1];
}

template <typename T>
inline const T& _2dArray<T>::operator()(size_t index0, size_t index1) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1);
#endif
    return m_storage[index0 + m_extents[0] * index1];
}

template <typename T>
inline size_t _2dArray<T>::extent(size_t dimension) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_dimension(dimension);
#endif
    return m_extents[dimension];
}

template <typename T>
inline void _2dArray<T>::set_size(size_t extent0, size_t extent1)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1);
#endif
    if (extent0 * extent1 == m_extents[0] * m_extents[1]) {
        m_extents[0] = extent0;
        m_extents[1] = extent1;
    }
    else {
        free_memory();
        init_memory(extent0, extent1);
    }
}

template <typename T>
inline typename _2dArray<T>::iterator _2dArray<T>::begin()
{
    return m_storage;
}

template <typename T>
inline typename _2dArray<T>::const_iterator _2dArray<T>::begin() const
{
    return m_storage;
}

template <typename T>
inline typename _2dArray<T>::iterator _2dArray<T>::end()
{
    return m_storage + m_extents[0] * m_extents[1];
}

template <typename T>
inline typename _2dArray<T>::const_iterator _2dArray<T>::end() const
{
    return m_storage + m_extents[0] * m_extents[1];
}

#ifdef FIBER_CHECK_ARRAY_BOUNDS
template <typename T>
inline void _2dArray<T>::check_dimension(size_t dimension) const
{
    if (dimension < 0 || 1 < dimension)
        throw std::invalid_argument("Invalid dimension");
}

template <typename T>
inline void _2dArray<T>::check_extents(size_t extent0, size_t extent1) const
{
    if (extent0 <= 0 || extent1 <= 0)
        throw std::length_error("Invalid extent");
}

template <typename T>
inline void _2dArray<T>::check_indices(size_t index0, size_t index1) const
{
    if (index0 < 0 || m_extents[0] <= index0 ||
        index1 < 0 || m_extents[1] <= index1)
        throw std::out_of_range("Invalid index");
}
#endif // FIBER_CHECK_ARRAY_BOUNDS

// _1dSliceOf2dArray

template <typename T>
inline _1dSliceOf2dArray<T>::_1dSliceOf2dArray(
        _2dArray<T>& array, size_t index1) :
    m_array(array), m_index1(index1)
{}

template <typename T>
inline _1dSliceOf2dArray<T>& _1dSliceOf2dArray<T>::self() {
    return *this;
}

template <typename T>
inline const T& _1dSliceOf2dArray<T>::operator()(size_t index0) const {
    return m_array(index0, m_index1);
}

template <typename T>
inline T& _1dSliceOf2dArray<T>::operator()(size_t index0) {
    return m_array(index0, m_index1);
}

template <typename T>
inline size_t _1dSliceOf2dArray<T>::extent(size_t dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _1dSliceOf2dArray<T>::check_dimension(size_t dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (dimension < 0 || 0 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _1dSliceOfConst2dArray

template <typename T>
inline _1dSliceOfConst2dArray<T>::_1dSliceOfConst2dArray(
        const _2dArray<T>& array, size_t index1) :
    m_array(array), m_index1(index1)
{}

template <typename T>
inline const T& _1dSliceOfConst2dArray<T>::operator()(size_t index0) const {
    return m_array(index0, m_index1);
}

template <typename T>
inline size_t _1dSliceOfConst2dArray<T>::extent(size_t dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _1dSliceOfConst2dArray<T>::check_dimension(size_t dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (dimension < 0 || 0 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

} // namespace Fiber
