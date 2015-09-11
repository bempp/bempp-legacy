from cython.operator cimport dereference as deref
from .discrete_boundary_operator cimport RealDiscreteBoundaryOperator
from .discrete_boundary_operator cimport ComplexDiscreteBoundaryOperator
from bempp_ext.utils.parameter_list cimport ParameterList
from bempp_ext.space.space cimport Space
from .assembler cimport c_LocalAssemblerForIntegralOperators
from .assembler cimport c_LocalAssemblerForLocalOperators
from .assembler cimport RealIntegralOperatorLocalAssembler
from .assembler cimport ComplexIntegralOperatorLocalAssembler
from .assembler cimport LocalOperatorLocalAssembler



cdef class RealElementaryIntegralOperator:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def make_local_assembler(self, ParameterList parameters):
        cdef RealIntegralOperatorLocalAssembler assembler = RealIntegralOperatorLocalAssembler()
        assembler.impl_ = deref(self.impl_).makeAssembler(deref(parameters.impl_))
        return assembler

    def assemble_weak_form(self, ParameterList parameters):
        cdef RealDiscreteBoundaryOperator op = RealDiscreteBoundaryOperator()
        op.impl_ = deref(self.impl_).assembleWeakForm(deref(parameters.impl_))
        return op

    property domain:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).domain()
            return space

    property range:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).range()
            return space

    property dual_to_range:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).dualToRange()
            return space

cdef class ComplexElementaryIntegralOperator:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def make_local_assembler(self, ParameterList parameters):
        cdef ComplexIntegralOperatorLocalAssembler assembler = ComplexIntegralOperatorLocalAssembler()
        assembler.impl_ = deref(self.impl_).makeAssembler(deref(parameters.impl_))
        return assembler

    def assemble_weak_form(self, ParameterList parameters):
        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()
        op.impl_ = deref(self.impl_).assembleWeakForm(deref(parameters.impl_))
        return op

    property domain:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).domain()
            return space

    property range:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).range()
            return space

    property dual_to_range:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).dualToRange()
            return space

cdef class ElementaryLocalOperator:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def make_local_assembler(self, ParameterList parameters):
        cdef LocalOperatorLocalAssembler assembler = LocalOperatorLocalAssembler()
        assembler.impl_ = deref(self.impl_).makeAssembler(deref(parameters.impl_))
        return assembler

    def assemble_weak_form(self, ParameterList parameters):
        cdef RealDiscreteBoundaryOperator op = RealDiscreteBoundaryOperator()
        op.impl_ = deref(self.impl_).assembleWeakForm(deref(parameters.impl_))
        return op
        
    property domain:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).domain()
            return space

    property range:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).range()
            return space

    property dual_to_range:
        def __get__(self):
            cdef Space space = Space()
            space.impl_ = deref(self.impl_).dualToRange()
            return space

