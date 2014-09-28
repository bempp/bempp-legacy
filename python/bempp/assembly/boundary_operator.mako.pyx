<%
from bempp_operators import dtypes, compatible_dtypes, bops
ifloop = lambda x: "if" if getattr(x, 'index', x) == 0 else "elif"
%>\
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from bempp.options cimport Options

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.values():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor

<%def name="type_loop(inner_text)" filter="trim">
% for i, (pybasis, cybasis) in enumerate(dtypes.items()):
%     for j, (pyresult, cyresult) in enumerate(dtypes.items()):
%         if pyresult in compatible_dtypes[pybasis]:
        ${ifloop(i + j)} self.basis_type == ${i} and self.result_type == ${j}:
            ${inner_text(cybasis, cyresult)}
%         endif
%     endfor
% endfor
        else:
            msg = "Unknown or incompatible basis and result types"
            raise RuntimeError(msg)
</%def>


cdef extern from "bempp/assembly/symmetry.hpp" namespace "Bempp":
    cdef enum Symmetry:
        NO_SYMMETRY
        SYMMETRIC
        HERMITIAN
        AUTO_SYMMETRY


cdef extern from "bempp/assembly/python.hpp":
    void inplace_boundary_operator[BASIS, RESULT](void* memory)
    void deconstructor[BASIS, RESULT](void* memory)



cdef class BoundaryOperator:
    """ Holds a reference to a boundary operator """

    def __cinit__(self, *args, **kwargs):
        """ Initializes memory """
        self.constructed = False
        cdef int n = max([
% for pybasis, cybasis in dtypes.items():
%     for pyresult, cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            sizeof(c_BoundaryOperator[${cybasis}, ${cyresult}]),
%         endif
%     endfor
% endfor
        ])
        memory = PyMem_Malloc(n)
        if not memory:
            self.memory = NULL
            raise MemoryError("Could not allocate memory")
        self.memory = memory

    def __init__(self, basis_type=None, result_type=None):
        from numpy import dtype

        # Simple way to check basis and result type are compatible
        Options(basis_type=basis_type, result_type=result_type)

<% type_numbers = {k: i for i, k in enumerate(dtypes)} %>
        self.basis_type = ${repr(type_numbers)}[str(basis_type)]
        self.result_type = ${repr(type_numbers)}[str(result_type)]

        # Call to in-place new:
<%
    def create_inplace(cybasis, cyresult):
        return "inplace_boundary_operator[%s, %s](self.memory)" % (
            cybasis, cyresult
        )
%>\
        ${type_loop(create_inplace)}

        self.constructed = True

    def __dealloc__(self):
        # indentation issues mean we can't put type loop within if
        # statement. So this function has three exit points:
        # 1. Memory was not allocated
        # 2. Memory was allocated but object was not constructed
        # 3. Memory was allocated and object was constructed
        if not self.memory:
            return

        cdef void * memory = self.memory
        constructed = self.constructed
        # disables object before deconstructing and deallocating
        self.memory = NULL
        self.constructed = False

        if not constructed:
            PyMem_Free(memory)
            return

        ${type_loop(lambda x, y:"deconstructor[%s, %s](memory)" % (x,y))}
        PyMem_Free(memory)

% for variable in ['basis', 'result']:
    property ${variable}_type:
        def __get__(self):
            from numpy import dtype
            return dtype(${repr(list(dtypes.keys()))}[self.${variable}_type])
% endfor
