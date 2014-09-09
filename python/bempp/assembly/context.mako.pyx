<%
from bempp_operators import bops, dtypes, compatible_dtypes
ifloop = lambda x: "if" if getattr(x, 'index', x) == 0 else "elif"
%>
from libcpp.string cimport string
from cython.operator cimport dereference as deref
from bempp.fiber.accuracy_options cimport c_AccuracyOptions
from bempp.assembly.boundary_operator cimport c_BoundaryOperator
from bempp.assembly.boundary_operator cimport BoundaryOperator
from bempp.options cimport AssemblyOptions, Options
from bempp.space.space cimport Space, c_Space
from bempp.utils cimport shared_ptr

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.itervalues():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor

cdef extern from "bempp/assembly/python.hpp":
% for opname, description in bops.iteritems():
%     if description['implementation'] == 'standard':
    c_BoundaryOperator[BASIS, RESULT] \
        ${description['c_creator']}[BASIS, RESULT](
            c_BoundaryOperator[BASIS, RESULT]& _output,
            const c_AccuracyOptions& accuracyOptions,
            const AssemblyOptions& assemblyOptions,
            const shared_ptr[c_Space[BASIS]]& domain,
            const shared_ptr[c_Space[BASIS]]& range,
            const shared_ptr[c_Space[BASIS]]& dualToRange,
            const string& label,
            int symmetry
        )
%     endif
% endfor

<%def name="create_operator(py_creator, c_creator, doc, **kwargs)">
    def ${py_creator}(self,
        Space domain, Space range, Space dual_to_range,
        str label="", int symmetry=0):
        """ ${doc} """

        # Check space types match this context
        if self.basis_type != range.dtype \
            or self.basis_type != dual_to_range.dtype \
            or self.basis_type != domain.dtype:
                raise TypeError("Incompatible spaces")

        # Creates C accuracy option
        cdef c_AccuracyOptions acc_ops
        self._accuracy.to_cpp(acc_ops)

        result = BoundaryOperator(basis_type=self.basis_type,
                result_type=self.result_type)

        # Loop over possible basis and result type combinations,
        # And call templated function to create the operator
% for i, (pybasis, cybasis) in enumerate(dtypes.iteritems()):
%    for j, (pyresult, cyresult) in enumerate(dtypes.iteritems()):
%       if pyresult in compatible_dtypes[pybasis]:
        ${ifloop(i + j)} self.basis_type == '${pybasis}' ${'\\'}
            and self.result_type == '${pyresult}':
                ${c_creator}[${cybasis}, ${cyresult}](
                    deref(
                        <c_BoundaryOperator[${cybasis}, ${cyresult}]*>
                        result.memory
                    ),
                    acc_ops, self.assembly,
                    domain.impl_${pybasis},
                    range.impl_${pybasis},
                    dual_to_range.impl_${pybasis},
                    label, symmetry
                )
%       endif
%    endfor
% endfor
        else:
            msg = "Unknown or incompatible basis and result types"
            raise RuntimeError(msg)

        return result
</%def>

cdef class Context(Options):
    cdef:
        ## simplified access to operators
        readonly object operators

    def __init__(self, **kwargs):
        from ..nested_attributes import NestedAttributes
        super(Context, self).__init__(**kwargs)
        self.operators = NestedAttributes({
% for key, description in bops.iteritems():
<%
    location = description['location']
    py_creator = description['py_creator']
%>\
            ${repr(location)}: self.${py_creator},
% endfor
        })
        """ All operators available to this context """

    def scalar_space(self, grid, *args, **kwargs):
        """ Creates a scalar space

            This is a convenience function for creating scalar spaces. It
            merely avoids refering to the `basis_type`. The input is the same
            as :py:func:`bempp.space.scalar_spaces`, but *without* the `dtypes`
            argument.
        """
        from .. import scalar_space
        return scalar_space(grid, self.basis_type, *args, **kwargs)


% for description in bops.itervalues():
%     if description['implementation'] == 'standard':
    ${create_operator(**description)}
%     endif
% endfor
