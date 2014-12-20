from libcpp cimport bool as cbool

cdef extern from "<armadillo>" namespace "arma":
    ctypedef size_t uword
    cdef cppclass Col[T]:
        Col(T*, uword, cbool copy_aux_mem, cbool strict) nogil
        Col() nogil
        T* memptr()
        T& value "operator()"(int i) # bounds checking
        T& at(int i) # No bounds checking
        int n_rows
        int n_cols

    cdef cppclass Mat[T]:
        Mat(T*, uword, uword, cbool copy_aux_mem, cbool strict) nogil
        Mat() nogil
        T* memptr()
