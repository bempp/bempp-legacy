from bempp.core.utils cimport unique_ptr
from bempp.core.utils cimport complex_double
from bempp.core.utils cimport Matrix
from libcpp.vector cimport vector

cdef extern from "bempp/fiber/local_assembler_for_integral_operators.hpp" namespace "Fiber":
    cdef cppclass c_LocalAssemblerForIntegralOperators "Fiber::LocalAssemblerForIntegralOperators"[T]:
        pass

cdef extern from "bempp/fiber/local_assembler_for_local_operators.hpp" namespace "Fiber":
    cdef cppclass c_LocalAssemblerForLocalOperators "Fiber::LocalAssemblerForLocalOperators"[T]:
        void evaluateLocalWeakForms(const vector[int]&, vector[Matrix[T]]&)


cdef class RealIntegralOperatorLocalAssembler:
    cdef unique_ptr[c_LocalAssemblerForIntegralOperators[double]] impl_

cdef class ComplexIntegralOperatorLocalAssembler:
    cdef unique_ptr[c_LocalAssemblerForIntegralOperators[complex_double]] impl_

cdef class LocalOperatorLocalAssembler:
    cdef unique_ptr[c_LocalAssemblerForLocalOperators[double]] impl_


