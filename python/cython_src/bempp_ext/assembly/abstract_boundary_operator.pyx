from cython.operator cimport dereference as deref
from bempp.utils.parameter_list cimport ParameterList

cdef class Assembler:
    pass

cdef class RealIntegralOperatorLocalAssembler(Assembler):
    pass

cdef class ComplexIntegralOperatorLocalAssembler(Assembler):
    pass

cdef class LocalOperatorAssembler(Assembler):
    pass

cdef class ElementaryBoundaryOperator:

    def make_assembler(ParameterList parameter_list):

        raise NotImplementedError("Method not implemented")

    def assemble_weak_form(ParameterList parameter_list):

        raise NotImplementedError("Method not implemented")

cdef class RealElementaryIntegralOperator(ElementaryBoundaryOperator):

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def make_assembler(self, ParameterList parameters):
        cdef RealIntegralOperatorLocalAssembler assembler = RealIntegralOperatorLocalAssembler()
        assembler.impl_ = deref(self.impl_).makeAssembler(deref(parameters.impl_))
        return assembler

cdef class ComplexElementaryIntegralOperator(ElementaryBoundaryOperator):

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def make_assembler(self, ParameterList parameters):
        cdef ComplexIntegralOperatorLocalAssembler assembler = ComplexIntegralOperatorLocalAssembler()
        assembler.impl_ = deref(self.impl_).makeAssembler(deref(parameters.impl_))
        return assembler

cdef class ElementaryLocalOperator(ElementaryBoundaryOperator):

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def make_assembler(self, ParameterList parameters):
        cdef LocalOperatorLocalAssembler assembler = LocalOperatorLocalAssembler()
        assembler.impl_ = deref(self.impl_).makeAssembler(deref(parameters.impl_))
        return assembler
        

