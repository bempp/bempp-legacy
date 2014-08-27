<%
from aca_options import properties, enums, enum_properties

def pythonic(name):
    name = str(name)
    to_lower = lambda x: x if x == x.lower() else "_" + x.lower() 
    return name[0].lower() + ''.join([to_lower(u) for u in name[1:]])

def enum_names(name):
    return name.lower()         \
        .replace("_", " ")      \
        .replace("mode", "")    \
        .replace("reaction", "")\
        .replace("assembly", "").rstrip().lstrip()

def ifloop(loop):
    return 'if' if loop.index == 0 else 'elif'

def allowables(name):
    return '|'.join([repr(enum_names(u)) for u in enums[name]])
%>

cdef extern from "<limits.h>":
    int UINT_MAX

cdef class Options:
    """ Handles all options for BEM++

        Parameters
        ----------
% for name, (vartype, default, doc) in properties.iteritems():

        ${pythonic(name)}: ${vartype}
            ${doc.rstrip().lstrip()}${'.' if doc.rstrip()[-1] != '.' else ''}
            Defaults to ${default}.
% endfor
% for name, (default, doc) in enum_properties.iteritems():

        ${pythonic(name)}: ${allowables(name)}
            ${doc.rstrip().lstrip()}${'.' if doc.rstrip()[-1] != '.' else ''}
            Defaults to ${default}.
% endfor

        do_opencl: bool
            Whether to use opencl. Defaults to False.

        max_threads: int or 'auto'
            Maximum number of threads used during assembly. Defaults to 'auto'.
    """
    def __init__(self, **kwargs):
<%
iterator = properties.keys() \
    + enum_properties.keys() \
    + ['do_opencl', 'max_threads']
%>
% for name in iterator:
        self.${pythonic(name)} = kwargs.get("${pythonic(name)}", None)
% endfor

% for name, (vartype, default, docstring) in properties.iteritems():
    property ${pythonic(name)}:
        """ ${docstring} """
        def __get__(self):
            return self.aca_options.${name}
        def __set__(self, value):
            self.aca_options.${name} = ${default} if value is None else value
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


% for name, values in enums.iteritems():
    <%
    cname = 'mode' if name == 'AcaAssemblyMode' \
            else 'reactionToUnsupportedMode'
    %>
    property ${pythonic(name)}:
        """ ${enum_properties[name][1]}

            Can take the following values: ${allowables(name)}
        """
        def __get__(self):
    % for value in values:
            ${ifloop(loop)} self.aca_options.${cname} == ${value}:
                return "${enum_names(value)}"
    % endfor
            else:
                raise RuntimeError("C value is invalid")
        def __set__(self, value):
            if value is None:
                self.${pythonic(name)} = "${enum_properties[name][0]}"
    % for value in values:
            elif value == "${enum_names(value)}":
                self.aca_options.${cname} = ${value}
    % endfor
            else:
                raise ValueError("Incorrect input %s" % value)
% endfor
