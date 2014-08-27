<%
from bempp_options import options

def ifloop(loop):
    return 'if' if loop.index == 0 else 'elif'

def getter_func(origin):
    condition = origin == 'AcaOptions'
    return 'aca_options.%s' if condition else 'assembly.%s()'
def setter_func(origin):
    condition = origin == 'AcaOptions'
    return 'aca_options.%s = %s' if condition else 'assembly.%s(%s)'

%>

<%def name="getsetter_impl(getter, default, setter, **kwargs)">
<%
    getter_string = getter_func(kwargs['c origin'])
    setter_string = setter_func(kwargs['c origin'])
%>
        def __get__(self):
            return self.${getter_string % getter}
        def __set__(self, value):
            if value is None:
                value = ${default}
            self.${setter_string % (setter, 'value')}
</%def>

<%def name="enums_impl(type, getter, setter, default, enums, **kwargs)">
<%
    getter_string = getter_func(kwargs['c origin'])
    setter_string = setter_func(kwargs['c origin'])
%>
        def __get__(self):
            cdef ${type} c_value = self.${getter_string % getter}
%   for python, cython in enums.iteritems():
            ${ifloop(loop)} c_value == ${cython}:
                return ${repr(python)}
%   endfor
            else:
                raise RuntimeError("C value is invalid")
        def __set__(self, value):
            if value is None:
                self.${setter_string % (setter, default)}
%   for python, cython in enums.iteritems():
            elif ${repr(python)} == value:
                self.${setter_string % (setter, cython)}
%   endfor
            else:
                raise ValueError("Incorrect input %s" % value)
</%def>

cdef extern from "<limits.h>":
    int UINT_MAX

cdef class Options:
    """ Handles all options for BEM++

        Parameters
        ----------
% for option, description in options.iteritems():
<%
    docstring = description['doc'].rstrip().lstrip()
    if docstring[-1] == '.': docstring = docstring[:-1]
%>
        ${description['pyname']}: ${description['doc_type']}
            ${docstring}. Defaults to ${description['doc_default']}.
% endfor

        do_opencl: bool
            Whether to use opencl. Defaults to False.

        max_threads: int or 'auto'
            Maximum number of threads used during assembly. Defaults to 'auto'.
    """
    def __init__(self, **kwargs):
% for name, description in options.iteritems():
<%
    pyname = description['pyname']
%>
        self.${pyname} = kwargs.get("${pyname}", None)
% endfor
        self.do_opencl = kwargs.get("do_opencl", None)
        self.max_threads = kwargs.get("max_threads", None)

% for option, description in options.iteritems():
    property ${description['pyname']}:
        """ ${description['doc']}

            Allowed input:
                ${description['doc_type']} or None (reverts to default)
        """
%  if description['implementation'] == 'getsetters':
${getsetter_impl(**description)}
%   elif description['implementation'] == 'enums':
${enums_impl(**description)}
%   endif
% endfor


    property do_opencl:
        def __get__(self):
            return self.parallelization.isOpenClEnabled()
        def __set__(self, value):
            cdef OpenClOptions opencl_ops
            if value is not None or value == True:
                self.parallelization.enableOpenCl(opencl_ops)
            else:
                self.parallelization.disableOpenCl()
    property max_threads:
        def __get__(self):
            if self.parallelization.maxThreadCount() < 0:
                return 'auto'
            return self.parallelization.maxThreadCount()
        def __set__(self, value):
            if value is None or value == 'auto' or value < 0:
                value = -1
            self.parallelization.setMaxThreadCount(value)
