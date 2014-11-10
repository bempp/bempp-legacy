<%
from bempp_options import options


def assign(argument):
    return "self.{0} = kwargs.pop('{0}', None)".format(argument)


def getter_func(origin):
    condition = origin == 'AcaOptions'
    return 'aca_options.%s' if condition else 'assembly.%s()'


def setter_func(origin):
    condition = origin == 'AcaOptions'
    return 'aca_options.%s = %s' if condition else 'assembly.%s(%s)'
%>


<%def name="getsetter_impl(getter, default, setter, **description)">
<%
    getter_string = getter_func(description['c origin'])
    setter_string = setter_func(description['c origin'])
%>
        def __get__(self):
            return self.${getter_string % getter}
        def __set__(self, value):
            if value is None:
                value = ${default}
            self.${setter_string % (setter, 'value')}
%   if description['c origin'] == "AcaOptions":
            # Aca options hidden inside assembly structure
            # So we must reset them all each time any one changes
            # Avoid this operation when not using Aca assembly
            if self.assembly.assemblyMode() == ASSEMBLY_MODE_ACA:
                self.assembly.switchToAcaMode(self.aca_options)
%   endif
</%def>


<%def name="enums_impl(type, getter, setter, default, enums, **description)">
<%
    getter_string = getter_func(description['c origin'])
    setter_string = setter_func(description['c origin'])
    ifloop = lambda x: 'if' if x.index == 0 else 'elif'
%>
        def __get__(self):
            cdef ${type} c_value = self.${getter_string % getter}
%   for python, cython in enums.items():
            ${ifloop(loop)} c_value == ${cython}:
                return ${repr(python)}
%   endfor
            else:
                raise RuntimeError("C value is invalid")
        def __set__(self, value):
            if value is None:
                self.${setter_string % (setter, default)}
%   for python, cython in enums.items():
            elif ${repr(python)} == value:
                self.${setter_string % (setter, cython)}
%   endfor
            else:
                raise ValueError("Incorrect input %s" % value)
</%def>


<%def name="automatic_impls(**description)" filter="trim">
%   if description['implementation'] != 'manual':
        ${property_and_doc(**description)}\
%       if description['implementation'] == 'getsetters':
            ${getsetter_impl(**description)}
%       elif description['implementation'] == 'enums':
            ${enums_impl(**description)}
%       endif
%   endif
</%def>


<%def name="property_and_doc(**description)" filter="trim">
    property ${description['pyname']}:
        """ ${description['doc']}

            Allowed input:
                ${description['doc_type']} or None (reverts to default)
        """
</%def>


cdef extern from "<limits.h>":
    int UINT_MAX

cdef class Options:
    """ Handles all options for BEM++

        Parameters
        ----------
% for option, description in options.items():
<%
    docstring = description['doc'].rstrip().lstrip()
    if docstring[-1] == '.': docstring = docstring[:-1]
%>
        ${description['pyname']}: ${description['doc_type']}
            ${docstring}. Defaults to ${description['doc_default']}.
% endfor

        basis_type: float or complex numpy.dtype
            Type used to represent the values of the (components of the) basis
            functions into which arguments of operators discretized with this
            strategy will be expanded. Defaults to 'float64'.

        result_type: float or complex numpy.dytpe
            Type used to represent the values of integrals.
            Defaults to 'complex128'.


        Note
        ----

        The precision (single or double) of the basis and result  must be
        strictly equal. The kind (real or imaginary) of the result must larger
        or equal to the kind of the basis. When modifying the size of one, the
        size of the other is adjusted automatically. When modifying the kind of
        one type, adjustements are made to other only if stricly necessary.
    """
    def __init__(self, **kwargs):
% for name, description in options.items():
        ${description['pyname'] | assign}
% endfor
        # Both basis and result type are given:
        # Check they are compatible
        basis = kwargs.pop('basis_type', None)
        result = kwargs.pop('result_type', None)
        if basis is not None and result is not None:
            self.result_type = result
            self.basis_type = basis
            if result != self.result_type:
                raise TypeError("Incompatible basis type")
        # Only basis type is set, so set it last
        elif basis is not None:
            self.result_type = None
            self.basis_type = basis
        # basis type is none, so set result type last
        else:
            self.basis_type = None
            self.result_type = result

        if len(kwargs) != 0:
            raise TypeError("Unexpected argument(s): %s" % kwargs.keys())

% for option, description in options.items():

    ${automatic_impls(**description)}
% endfor

    ${property_and_doc(**options['max_threads'])}
        def __get__(self):
            return 'auto' if self._max_threads < 0 else self._max_threads
        def __set__(self, value):
            cdef cbool cond = value is None or value == 'auto' or value < 0
            self._max_threads = -1 if cond else value
            self.assembly.setMaxThreadCount(self._max_threads)

    ${property_and_doc(**options['assembly_mode'])}
        def __get__(self):
            cdef AssemblyMode mode = self.assembly.assemblyMode()
            if mode == ASSEMBLY_MODE_ACA:
                return "aca"
            elif mode == ASSEMBLY_MODE_DENSE:
                return "dense"
            else:
                msg = "Internal bug: unknown assembly mode %s" % mode
                raise RuntimeError(msg)
        def __set__(self, value):
            if value is None:
                value = "dense"
            if value == "dense":
                self.assembly.switchToDenseMode()
            elif value == "aca":
                self.assembly.switchToAcaMode(self.aca_options)

    ${property_and_doc(**options['accuracy'])}
        def __get__(self):
            return self._accuracy
        def __set__(self, value):
            self._accuracy = AccuracyOptions(value) if value is not None \
                else AccuracyOptions()

    property basis_type:
        def __get__(self):
            return self._basis_type
        def __set__(self, value):
            from numpy import dtype
            value = str(value) if value is not None else 'float64'
            if value not in ['float32', 'float64', 'complex64', 'complex128']:
                raise ValueError("Incorrect basis type (%s)" % value)
            self._basis_type = dtype(value)

            # Promote result type accordingly
            cdef:
                size = numpy_size(value)
                kind = numpy_kind(value) and numpy_kind(str(self._result_type))
            self._result_type = dtype(get_numpy_type(size, kind))


    property result_type:
        def __get__(self):
            return self._result_type
        def __set__(self, value):
            from numpy import dtype
            value = str(value) if value is not None else 'float64'
            if value not in ['float32', 'float64', 'complex64', 'complex128']:
                raise ValueError("Incorrect basis type (%s)" % value)
            self._result_type = dtype(value)

            # Promote basis type accordingly
            cdef:
                size = numpy_size(value)
                kind = numpy_kind(value) or numpy_kind(str(self._basis_type))
            self._basis_type = dtype(get_numpy_type(size, kind))

    def __str__(self):
        return str({
%   for name, description in options.items():
            '${description['pyname']}': self.${description['pyname']},
%   endfor
        })

    def __repr__(self):
        # Figures out defaults so we print the minimum amount of necessary
        # information
        defaults = {
%   for name, description in options.items():
            '${description['pyname']}': ${description['doc_default']},
%   endfor
        }
        kwargs = {
            k: getattr(self, k) for k, default in defaults.items()
                if default != getattr(self, k)
        }
        return "bempp.Options(**%s)" % kwargs if len(kwargs) \
                else "bempp.Options()"

cdef int numpy_size(str x):
    return 4 if x == 'float32' or x == 'complex64' else 8

cdef cbool numpy_kind(str x):
    return x == 'float32' or x == 'float64'

cdef str get_numpy_type(int size, cbool is_real):
    if size == 4:
        return 'float32' if is_real else 'complex64'
    return 'float64' if is_real else 'complex128'
