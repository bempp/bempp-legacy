#include <stdexcept>

namespace Fiber
{

template <typename T>
inline Array3D<T>::Array3D()
{
    m_extents[0] = 0;
    m_extents[1] = 0;
    m_extents[2] = 0;
    m_storage = 0;
    m_owns = false;
}

template <typename T>
inline Array3D<T>::Array3D(int extent0, int extent1, int extent2)
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
inline Array3D<T>::Array3D(int extent0, int extent1, int extent2, T* data)
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
inline Array3D<T>::~Array3D()
{
    if (m_owns && m_storage)
        delete[] m_storage;
}

template <typename T>
inline T& Array3D<T>::operator()(int index0, int index1, int index2)
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
inline const T& Array3D<T>::operator()(int index0, int index1, int index2) const
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
inline int Array3D<T>::extent(int dimension) const
{
#ifdef FIBER_CHECK_ARRAY_BOUNDS
    check_dimension(dimension);
#endif
    return m_extents[dimension];
}

template <typename T>
inline void Array3D<T>::set_size(int extent0, int extent1, int extent2)
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
inline typename Array3D<T>::iterator Array3D<T>::begin()
{
    return m_storage;
}

template <typename T>
inline typename Array3D<T>::const_iterator Array3D<T>::begin() const
{
    return m_storage;
}

template <typename T>
inline typename Array3D<T>::iterator Array3D<T>::end()
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2];
}

template <typename T>
inline typename Array3D<T>::const_iterator Array3D<T>::end() const
{
    return m_storage + m_extents[0] * m_extents[1] * m_extents[2];
}

#ifdef FIBER_CHECK_ARRAY_BOUNDS
template <typename T>
inline void Array3D<T>::check_dimension(int dimension) const
{
    if (dimension < 0 || 2 < dimension)
        throw std::invalid_argument("Invalid dimension");
}

template <typename T>
inline void Array3D<T>::check_extents(int extent0, int extent1, int extent2) const
{
    if (extent0 <= 0 || extent1 <= 0 || extent2 <= 0)
        throw std::length_error("Invalid extent");
}

template <typename T>
inline void Array3D<T>::check_indices(int index0, int index1, int index2) const
{
    if (index0 < 0 || m_extents[0] <= index0 ||
        index1 < 0 || m_extents[1] <= index1 ||
        index2 < 0 || m_extents[2] <= index2)
        throw std::out_of_range("Invalid index");
}
#endif // FIBER_CHECK_ARRAY_BOUNDS

} // namespace Fiber
