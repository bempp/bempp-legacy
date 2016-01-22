

cdef extern from "bempp/common/multidimensional_arrays.hpp" namespace "Fiber":
    cdef cppclass _4dArray "Fiber::_4dArray<double>":
        _4dArray()
        size_t extent(size_t dimension) const
        const double* begin() const

    cdef cppclass _3dArray "Fiber::_3dArray<double>":
        _3dArray()
        size_t extent(size_t dimension) const
        const double* begin() const

cdef object _3d_array_to_numpy(_3dArray)
cdef object _4d_array_to_numpy(_4dArray)
