from cython.operator cimport dereference as deref
from .discrete_boundary_operator cimport RealDiscreteBoundaryOperator
from .discrete_boundary_operator cimport ComplexDiscreteBoundaryOperator
from bempp.core.utils.parameter_list cimport ParameterList
from bempp.core.utils import _convert_to_bytes
from bempp.core.space.space cimport Space
from .assembler cimport c_LocalAssemblerForIntegralOperators
from .assembler cimport c_LocalAssemblerForLocalOperators
from .assembler cimport RealIntegralOperatorLocalAssembler
from .assembler cimport ComplexIntegralOperatorLocalAssembler
from .assembler cimport LocalOperatorLocalAssembler
from bempp.core.utils.enum_types cimport SymmetryMode, symmetry_mode
from bempp.core.fiber.shape_transformation_functors cimport ShapeTransformationFunctorContainer
from bempp.core.fiber.local_integrand_functors cimport LocalIntegrandFunctorContainer


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

def abstract_local_operator_from_functors_ext(Space domain, Space range_, Space dual_to_range,
        ShapeTransformationFunctorContainer test_functor,
        ShapeTransformationFunctorContainer trial_functor,
        LocalIntegrandFunctorContainer integrand_functor,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()

    op.impl_.assign(
            c_abstract_local_operator_from_functors(
                domain.impl_,
                range_.impl_,
                dual_to_range.impl_,
                _convert_to_bytes(label),
                symmetry_mode(_convert_to_bytes(symmetry)),
                deref(test_functor.impl_),
                deref(trial_functor.impl_),
                deref(integrand_functor.impl_)))
    return op



