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

#include "../common/common.hpp"

#include <stdexcept>

namespace Fiber
{

template <typename T>
inline Array4d<T>::Array4d()
{
    m_extents[0] = 0;
    m_extents[1] = 0;
    m_extents[2] = 0;
    m_extents[3] = 0;
    m_storage = 0;
    m_owns = false;
}

template <typename T>
inline Array4d<T>::Array4d(size_t extent0, size_t extent1, size_t extent2, size_t extent3)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2, extent3);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_extents[3] = extent3;
    m_storage = new T[extent0 * extent1 * extent2 * extent3];
    m_owns = true;
}

template <typename T>
inline Array4d<T>::Array4d(size_t extent0, size_t extent1, size_t extent2, size_t extent3, T* data)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2, extent3);
#endif
    m_extents[0] = extent0;
    m_extents[1] = extent1;
    m_extents[2] = extent2;
    m_extents[3] = extent3;
    m_storage = data;
    m_owns = false;
}

template <typename T>
inline Array4d<T>::~Array4d()
{
    if (m_owns && m_storage)
        delete[] m_storage;
}

template <typename T>
inline T& Array4d<T>::operator()(size_t index0, size_t index1, size_t index2, size_t index3)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1, index2, index3);
#endif
    return m_storage[
            index0 +
            m_extents[0] * index1 +
            m_extents[0] * m_extents[1] * index2 +
            m_extents[0] * m_extents[1] * m_extents[2] * index3];
}

template <typename T>
inline const T& Array4d<T>::operator()(size_t index0, size_t index1, size_t index2, size_t index3) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_indices(index0, index1, index2, index3);
#endif
    return m_storage[
            index0 +
            m_extents[0] * index1 +
            m_extents[0] * m_extents[1] * index2 +
            m_extents[0] * m_extents[1] * m_extents[2] * index3];
}

template <typename T>
inline size_t Array4d<T>::extent(size_t dimension) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_dimension(dimension);
#endif
    return m_extents[dimension];
}

template <typename T>
inline void Array4d<T>::set_size(size_t extent0, size_t extent1, size_t extent2, size_t extent3)
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_extents(extent0, extent1, extent2, extent3);
#endif
    if (extent0 * extent1 * extent2 * extent3 ==
            m_extents[0] * m_extents[1] * m_extents[2] * m_extents[3]) {
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
        m_extents[3] = extent3;
    }
    else {
        if (m_owns && m_storage) {
            delete[] m_storage;
            m_storage = 0;
        }
        m_extents[0] = extent0;
        m_extents[1] = extent1;
        m_extents[2] = extent2;
        m_extents[3] = extent3;
        m_storage = new T[extent0 * extent1 * extent2 * extent3];
        m_owns = true;
    }
}

template <typename T>
inline typename Array4d<T>::iterator Array4d<T>::begin()
{
    return m_storage;
}

template <typename T>
inline typename Array4d<T>::const_iterator Array4d<T>::begin() const
{
    return m_storage;
}

template <typename T>
inline typename Array4d<T>::iterator Array4d<T>::end()
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2] * m_extents[3];
}

template <typename T>
inline typename Array4d<T>::const_iterator Array4d<T>::end() const
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2] * m_extents[3];
}

#ifdef FIBER_CHECK_ARRAY_BOUNDS
template <typename T>
inline void Array4d<T>::check_dimension(size_t dimension) const
{
    if (dimension < 0 || 3 < dimension)
        throw std::invalid_argument("Invalid dimension");
}

template <typename T>
inline void Array4d<T>::check_extents(size_t extent0, size_t extent1, size_t extent2, size_t extent3) const
{
    if (extent0 <= 0 || extent1 <= 0 || extent2 <= 0 || extent3 <= 0)
        throw std::length_error("Invalid extent");
}

template <typename T>
inline void Array4d<T>::check_indices(size_t index0, size_t index1, size_t index2, size_t index3) const
{
    if (index0 < 0 || m_extents[0] <= index0 ||
        index1 < 0 || m_extents[1] <= index1 ||
        index2 < 0 || m_extents[2] <= index2 ||
        index3 < 0 || m_extents[3] <= index3)
        throw std::out_of_range("Invalid index");
}
#endif // FIBER_CHECK_ARRAY_BOUNDS

} // namespace Fiber
