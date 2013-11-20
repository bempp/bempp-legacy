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
    init_empty();
}

template <typename T>
inline _3dArray<T>::_3dArray(size_t extent0, size_t extent1, size_t extent2)
{
    init_memory(extent0, extent1, extent2);
}

template <typename T>
inline _3dArray<T>::_3dArray(size_t extent0, size_t extent1, size_t extent2, T* data, bool strict)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_storage = data;
    m_owns = false;
    m_strict = strict;
}

template <typename T>
inline _3dArray<T>::_3dArray(const _3dArray& other)
{
    if (!other.is_empty()) {
        init_memory(other.m_extents[0], other.m_extents[1], other.m_extents[2]);
        std::copy(other.begin(), other.end(), m_storage);
    } else
        init_empty();
}

template <typename T>
inline _3dArray<T>& _3dArray<T>::operator=(const _3dArray& rhs)
{
    if (&rhs != this) {
        set_size(rhs.m_extents[0], rhs.m_extents[1], rhs.m_extents[2]);
        std::copy(rhs.begin(), rhs.end(), m_storage);
    }
    return *this;
}

template <typename T>
inline _3dArray<T>::~_3dArray()
{
    free_memory();
}

template <typename T>
inline void _3dArray<T>::init_memory(size_t extent0, size_t extent1, size_t extent2)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2);
#endif
    m_storage = new T[extent0 * extent1 * extent2];
    m_owns = true;
    m_strict = false;
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
}

template <typename T>
inline void _3dArray<T>::init_empty()
{
    m_extents[0] = 0;
    m_extents[1] = 0;
    m_extents[2] = 0;
    m_storage = 0;
    m_owns = false;
    m_strict = false;
}

template <typename T>
inline void _3dArray<T>::free_memory()
{
    if (m_owns && m_storage)
        delete[] m_storage;
    m_owns = false;
    m_storage = 0;
}

template <typename T>
inline T& _3dArray<T>::operator()(size_t index0, size_t index1, size_t index2)
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
inline const T& _3dArray<T>::operator()(size_t index0, size_t index1, size_t index2) const
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
inline size_t _3dArray<T>::extent(size_t dimension) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_dimension(dimension);
#endif
    return m_extents[dimension];
}

template <typename T>
inline _3dArray<T>& _3dArray<T>::operator+=(const _3dArray<T>& other){
    if ((this->extent(0)!= other.extent(0))||
            (this->extent(1)!= other.extent(1))||
            (this->extent(2)!= other.extent(2)))
        std::runtime_error("_3dArray<T> operator+=: Array sizes don't agree.");
    for (size_t i=0;i<this->extent(2);++i)
        for (size_t j=0;j<this->extent(1);++j)
            for (size_t k=0;k<this->extent(0);++k)
                (*this)(k,j,i)+=other(k,j,i);
    return *this;

}

template <typename T>
inline _3dArray<T>& _3dArray<T>::operator*=(const T& other) {
    for (size_t i=0;i<this->extent(2);++i)
        for (size_t j=0;j<this->extent(1);++j)
            for (size_t k=0;k<this->extent(0);++k)
                (*this)(k,j,i)*=other;
    return *this;
}



template <typename T>
inline void _3dArray<T>::set_size(size_t extent0, size_t extent1, size_t extent2)
{
    if (extent0 * extent1 * extent2 ==
            m_extents[0] * m_extents[1] * m_extents[2]) {
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
    }
    else {
        if (m_strict)
            throw std::runtime_error("_3dArray::set_size(): Changing the total "
                                     "number of elements stored in an array "
                                     "created in the strict mode is not allowed");
        if (m_owns) {
            delete[] m_storage;
            m_storage = 0;
        }
        if (extent0 * extent1 * extent2 != 0)
            init_memory(extent0, extent1, extent2);
        else
            init_empty();
    }
}

template <typename T>
inline bool _3dArray<T>::is_empty() const
{
    return m_extents[0] * m_extents[1] * m_extents[2] == 0;
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
inline void _3dArray<T>::check_dimension(size_t dimension) const
{
    if (2 < dimension)
        throw std::invalid_argument("Invalid dimension");
}

template <typename T>
inline void _3dArray<T>::check_extents(size_t extent0, size_t extent1, size_t extent2) const
{
}

template <typename T>
inline void _3dArray<T>::check_indices(size_t index0, size_t index1, size_t index2) const
{
    if (m_extents[0] <= index0 ||
        m_extents[1] <= index1 ||
        m_extents[2] <= index2)
        throw std::out_of_range("Invalid index");
}
#endif // FIBER_CHECK_ARRAY_BOUNDS

// _2dSliceOf3dArray

template <typename T>
inline _2dSliceOf3dArray<T>::_2dSliceOf3dArray(_3dArray<T>& array, size_t index2) :
    m_array(array), m_index2(index2)
{}

template <typename T>
inline _2dSliceOf3dArray<T>& _2dSliceOf3dArray<T>::self() {
    return *this;
}

template <typename T>
inline const T& _2dSliceOf3dArray<T>::operator()(size_t index0, size_t index1) const {
    return m_array(index0, index1, m_index2);
}

template <typename T>
inline T& _2dSliceOf3dArray<T>::operator()(size_t index0, size_t index1) {
    return m_array(index0, index1, m_index2);
}

template <typename T>
inline size_t _2dSliceOf3dArray<T>::extent(size_t dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _2dSliceOf3dArray<T>::check_dimension(size_t dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (1 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _2dSliceOfConst3dArray

template <typename T>
inline _2dSliceOfConst3dArray<T>::_2dSliceOfConst3dArray(
        const _3dArray<T>& array, size_t index2) :
    m_array(array), m_index2(index2)
{}

template <typename T>
inline const T& _2dSliceOfConst3dArray<T>::operator()(size_t index0, size_t index1) const {
    return m_array(index0, index1, m_index2);
}

template <typename T>
inline size_t _2dSliceOfConst3dArray<T>::extent(size_t dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _2dSliceOfConst3dArray<T>::check_dimension(size_t dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (1 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _1dSliceOf3dArray

template <typename T>
inline _1dSliceOf3dArray<T>::_1dSliceOf3dArray(
        _3dArray<T>& array, size_t index1, size_t index2) :
    m_array(array), m_index1(index1), m_index2(index2)
{}

template <typename T>
inline _1dSliceOf3dArray<T>& _1dSliceOf3dArray<T>::self() {
    return *this;
}

template <typename T>
inline const T& _1dSliceOf3dArray<T>::operator()(size_t index0) const {
    return m_array(index0, m_index1, m_index2);
}

template <typename T>
inline T& _1dSliceOf3dArray<T>::operator()(size_t index0) {
    return m_array(index0, m_index1, m_index2);
}

template <typename T>
inline size_t _1dSliceOf3dArray<T>::extent(size_t dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _1dSliceOf3dArray<T>::check_dimension(size_t dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (0 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

// _1dSliceOfConst3dArray

template <typename T>
inline _1dSliceOfConst3dArray<T>::_1dSliceOfConst3dArray(
        const _3dArray<T>& array, size_t index1, size_t index2) :
    m_array(array), m_index1(index1), m_index2(index2)
{}

template <typename T>
inline const T& _1dSliceOfConst3dArray<T>::operator()(size_t index0) const {
    return m_array(index0, m_index1, m_index2);
}

template <typename T>
inline size_t _1dSliceOfConst3dArray<T>::extent(size_t dimension) const {
    check_dimension(dimension);
    return m_array.extent(dimension);
}

template <typename T>
inline void _1dSliceOfConst3dArray<T>::check_dimension(size_t dimension) const {
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    if (0 < dimension)
        throw std::invalid_argument("Invalid dimension");
#endif
}

} // namespace Fiber
