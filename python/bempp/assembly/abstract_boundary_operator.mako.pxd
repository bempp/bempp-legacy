from bempp.utils.unique_ptr cimport unique_ptr
from bempp.utils.shared_ptr cimport shared_ptr
from bempp.utils cimport complex_double
from bempp.utils.parameter_list cimport c_ParameterList
from bempp.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator

cdef extern from "bempp/fiber/local_assembler_for_integral_operators.hpp" namespace "Fiber":
    cdef cppclass c_LocalAssemblerForIntegralOperators "Fiber::LocalAssemblerForIntegralOperators"[T]:
        pass

cdef extern from "bempp/fiber/local_assembler_for_local_operators.hpp" namespace "Fiber":
    cdef cppclass c_LocalAssemblerForLocalOperators "Fiber::LocalAssemblerForLocalOperators"[T]:
        pass


cdef extern from "bempp/assembly/elementary_integral_operator.hpp" namespace "Bempp":
    cdef cppclass c_RealElementaryIntegralOperator "Bempp::ElementaryIntegralOperator<double,double,double>":
        unique_ptr[c_LocalAssemblerForIntegralOperators[double]] makeAssembler(const c_ParameterList&)
        shared_ptr[c_DiscreteBoundaryOperator[double]] assembleWeakForm(const c_ParameterList&)

    cdef cppclass c_ComplexElementaryIntegralOperator "Bempp::ElementaryIntegralOperator<double,std::complex<double>,std::complex<double> >":
        unique_ptr[c_LocalAssemblerForIntegralOperators[complex_double]] makeAssembler(const c_ParameterList&)
        shared_ptr[c_DiscreteBoundaryOperator[complex_double]] assembleWeakForm(const c_ParameterList&)

cdef extern from "bempp/assembly/elementary_local_operator.hpp" namespace "Bempp":
    cdef cppclass c_ElementaryLocalOperator "Bempp::ElementaryLocalOperator<double, double>":
        unique_ptr[c_LocalAssemblerForLocalOperators[double]] makeAssembler(const c_ParameterList&)
        shared_ptr[c_DiscreteBoundaryOperator[double]] assembleWeakForm(const c_ParameterList&)


cdef class Assembler:
    pass

cdef class RealIntegralOperatorLocalAssembler(Assembler):
    cdef unique_ptr[c_LocalAssemblerForIntegralOperators[double]] impl_

cdef class ComplexIntegralOperatorLocalAssembler(Assembler):
    cdef unique_ptr[c_LocalAssemblerForIntegralOperators[complex_double]] impl_

cdef class LocalOperatorLocalAssembler(Assembler):
    cdef unique_ptr[c_LocalAssemblerForLocalOperators[double]] impl_

cdef class ElementaryBoundaryOperator:
    pass

cdef class RealElementaryIntegralOperator(ElementaryBoundaryOperator):
    cdef shared_ptr[const c_RealElementaryIntegralOperator] impl_

cdef class ComplexElementaryIntegralOperator(ElementaryBoundaryOperator):
    cdef shared_ptr[const c_ComplexElementaryIntegralOperator] impl_

cdef class ElementaryLocalOperator(ElementaryBoundaryOperator):
    cdef shared_ptr[const c_ElementaryLocalOperator] impl_


