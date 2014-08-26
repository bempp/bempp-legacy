<%
from aca_options import properties, enums

def pythonic(name):
    name = str(name)
    to_lower = lambda x: x if x == x.lower() else "_" + x.lower() 
    return name[0].lower() + ''.join([to_lower(u) for u in name[1:]])
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
    """
    def __init__(self, **kwargs):
% for name, (vartype, default, docstring) in properties.iteritems():
        self.${name} = kwargs.get("${pythonic(name)}", ${default})
% endfor

% for name, (vartype, default, docstring) in properties.iteritems():
    property ${pythonic(name)}:
        """ ${docstring} """
        def __get__(self):
            return self.${name}
        def __set__(self, value):
            self.${name} = ${default} if value is None else value
% endfor


<%
def enum_names(name):
    return name.lower()         \
        .replace("_", " ")      \
        .replace("mode", "")    \
        .replace("reaction", "")\
        .replace("assembly", "").rstrip().lstrip()
%>
% for name, values in enums.iteritems():
    property ${pythonic(name)}:
        def __get__(self):
            return {
    % for value in values:
                ${value}: "${enum_names(value)}",
    % endfor
            }[self.${'_' + pythonic(name)}]
        def __set__(self, value):
            self.${'_' + pythonic(name)} = {
    % for value in values:
                "${enum_names(value)}": ${value},
    % endfor
                None: GLOBAL_ASSEMBLY
            }[value]
% endfor
