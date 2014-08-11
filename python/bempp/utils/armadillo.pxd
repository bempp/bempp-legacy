from libcpp cimport bool as cbool

cdef extern from "<armadillo>" namespace "arma":
    ctypedef size_t uword
    cdef cppclass Col[T]:
        Col(T*, uword, cbool copy_aux_mem, cbool strict) nogil
        Col() nogil
        T* memptr()

    cdef cppclass Mat[T]:
        Mat(T*, uword, uword, cbool copy_aux_mem, cbool strict) nogil
        Mat() nogil
        T* memptr()
