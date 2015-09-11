from bempp_ext.utils cimport unique_ptr
from bempp_ext.utils cimport complex_double

cdef extern from "bempp/fiber/local_assembler_for_integral_operators.hpp" namespace "Fiber":
    cdef cppclass c_LocalAssemblerForIntegralOperators "Fiber::LocalAssemblerForIntegralOperators"[T]:
        pass

cdef extern from "bempp/fiber/local_assembler_for_local_operators.hpp" namespace "Fiber":
    cdef cppclass c_LocalAssemblerForLocalOperators "Fiber::LocalAssemblerForLocalOperators"[T]:
        pass

cdef class LocalAssembler:
    pass

cdef class RealIntegralOperatorLocalAssembler(LocalAssembler):
    cdef unique_ptr[c_LocalAssemblerForIntegralOperators[double]] impl_

cdef class ComplexIntegralOperatorLocalAssembler(LocalAssembler):
    cdef unique_ptr[c_LocalAssemblerForIntegralOperators[complex_double]] impl_

cdef class LocalOperatorLocalAssembler(LocalAssembler):
    cdef unique_ptr[c_LocalAssemblerForLocalOperators[double]] impl_


