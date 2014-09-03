<%
from space import dtypes, compatible_dtypes
ifloop = lambda x: "if" if getattr(x, 'index', x) == 0 else "elif"
%>
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from bempp.options cimport Options

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.itervalues():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor

cdef extern from "bempp/assembly/symmetry.hpp" namespace "Bempp":
    cdef enum Symmetry:
        NO_SYMMETRY
        SYMMETRIC
        HERMITIAN
        AUTO_SYMMETRY

cdef extern from "bempp/assembly/boundary_operator.hpp":
    cdef cppclass c_BoundaryOperator "Bempp::BoundaryOperator" [BASIS, RESULT]:
        c_BoundaryOperator()

cdef extern from "bempp/assembly/python.hpp":
    void inplace_boundary_operator[BASIS, RESULT](void* memory)


cdef class BoundaryOperator:
    """ Holds a reference to a boundary operator """

    def __cinit__(self, *args, **kwargs):
        """ Initializes memory """
        cdef int n = max([
% for pybasis, cybasis in dtypes.iteritems():
%     for pyresult, cyresult in dtypes.iteritems():
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

        if basis_type not in ${repr(dtypes.keys() + [None])}:
            raise ValueError("Incorrect basis type")
        if result_type not in ${repr(dtypes.keys() + [None])}:
            raise ValueError("Incorrect result type")

        if basis_type is None:
            basis_type = Options(result_type=result_type).basis_type
        elif result_type is None:
            result_type = Options(basis_type=basis_type).result_type
        elif basis_type is not None and result_type is not None:
            options = Options(basis_type=basis_type, result_type=result_type)
            if basis_type != options.basis_type \
                or result_type != options.result_type:
                raise ValueError("Incompatible basis and result type")

<% type_numbers = {k: i for i, k in enumerate(dtypes)} %>
        self.basis_type = ${repr(type_numbers)}[basis_type]
        self.result_type = ${repr(type_numbers)}[result_type]

        # Call to in-place new:
% for i, (pybasis, cybasis) in enumerate(dtypes.iteritems()):
%     for j, (pyresult, cyresult) in enumerate(dtypes.iteritems()):
%         if pyresult in compatible_dtypes[pybasis]:
        ${ifloop(i + j)} self.basis_type != ${i} and self.result_type == ${j}:
            inplace_boundary_operator[${cybasis}, ${cyresult}](self.memory)
%         endif
%     endfor
% endfor

    def __dealloc__(self):
        cdef void * memory = self.memory
        self.memory = NULL
        PyMem_Free(memory)

% for variable in ['basis', 'result']:
    property ${variable}_type:
        def __get__(self):
            from numpy import dtype
            return dtype(${repr(dtypes.keys())}[self.${variable}_type])
% endfor
