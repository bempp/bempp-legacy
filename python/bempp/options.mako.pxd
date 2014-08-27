from libcpp cimport bool as cbool
<%
from aca_options import properties, enums

def transform(type):
    return {'bool': 'cbool'}.get(type, type)

def pythonic(name):
    name = str(name)
    to_lower = lambda x: x if x == x.lower() else "_" + x.lower() 
    return name[0].lower() + ''.join([to_lower(u) for u in name[1:]])
%>

# Declare enums from various files
cdef extern from "bempp/assembly/aca_options.hpp":
% for name, values in enums.iteritems():
    cdef enum ${name} "Bempp::AcaOptions::${name}":
    % for value in values:
        ${value} "Bempp::AcaOptions::${value}"
    % endfor
% endfor

# Declare c structure used in conversions
cdef extern from "bempp/assembly/aca_options.hpp" namespace "Bempp":
    cdef cppclass AcaOptions:
        AcaOptions()
% for name, (ptype, default, doc) in properties.iteritems():
        ${transform(ptype)} ${name}
% endfor
        ReactionToUnsupportedMode reactionToUnsupportedMode
        AcaAssemblyMode mode


cdef class Options:
    cdef:
% for name, (vartype, default, docstring) in properties.iteritems():
        ${transform(vartype)} ${name}
% endfor
% for name in enums.iterkeys():
        ${name} ${'_' + pythonic(name)}
% endfor
    cdef to_aca_options(self, AcaOptions *c_options)
