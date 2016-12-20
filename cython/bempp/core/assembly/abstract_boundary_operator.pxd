cimport numpy
from bempp.core.utils cimport unique_ptr
from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport complex_double
from bempp.core.utils cimport c_ParameterList
from bempp.core.space cimport c_Space
from .discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from .assembler cimport c_LocalAssemblerForIntegralOperators
from .assembler cimport c_LocalAssemblerForLocalOperators


cdef extern from "bempp/assembly/elementary_integral_operator.hpp" namespace "Bempp":
    cdef cppclass c_RealElementaryIntegralOperator "Bempp::ElementaryIntegralOperator<double,double,double>":
        unique_ptr[c_LocalAssemblerForIntegralOperators[double]] makeAssembler(const c_ParameterList&)
        shared_ptr[c_DiscreteBoundaryOperator[double]] assembleWeakForm(const c_ParameterList&)
        shared_ptr[const c_Space[double]] domain()
        shared_ptr[const c_Space[double]] range()
        shared_ptr[const c_Space[double]] dualToRange()

    cdef cppclass c_ComplexElementaryIntegralOperator "Bempp::ElementaryIntegralOperator<double,std::complex<double>,std::complex<double> >":
        unique_ptr[c_LocalAssemblerForIntegralOperators[complex_double]] makeAssembler(const c_ParameterList&)
        shared_ptr[c_DiscreteBoundaryOperator[complex_double]] assembleWeakForm(const c_ParameterList&)
        shared_ptr[const c_Space[double]] domain()
        shared_ptr[const c_Space[double]] range()
        shared_ptr[const c_Space[double]] dualToRange()

cdef extern from "bempp/assembly/elementary_local_operator.hpp" namespace "Bempp":
    cdef cppclass c_ElementaryLocalOperator "Bempp::ElementaryLocalOperator<double, double>":
        unique_ptr[c_LocalAssemblerForLocalOperators[double]] makeAssembler(const c_ParameterList&)
        shared_ptr[c_DiscreteBoundaryOperator[double]] assembleWeakForm(const c_ParameterList&)
        shared_ptr[const c_Space[double]] domain()
        shared_ptr[const c_Space[double]] range()
        shared_ptr[const c_Space[double]] dualToRange()



cdef class RealElementaryIntegralOperator:
    cdef shared_ptr[const c_RealElementaryIntegralOperator] impl_

cdef class ComplexElementaryIntegralOperator:
    cdef shared_ptr[const c_ComplexElementaryIntegralOperator] impl_

cdef class ElementaryLocalOperator:
    cdef shared_ptr[const c_ElementaryLocalOperator] impl_


