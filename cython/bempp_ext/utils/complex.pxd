cdef extern from "<complex>" namespace "std":
    cdef cppclass cpp_complex "std::complex"[T]:
        cpp_complex()
        cpp_complex(T alpha,T beta)

ctypedef cpp_complex[double] complex_double
