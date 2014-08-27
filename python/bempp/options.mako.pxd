from libcpp cimport bool as cbool
<%
from aca_options import properties, enums

def transform(type):
    return {'bool': 'cbool'}.get(type, type)
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

cdef extern from "bempp/fiber/opencl_options.hpp" namespace "Fiber":
    cdef cppclass OpenClOptions:
        OpenClOptions()

cdef extern from "bempp/fiber/parallelization_options.hpp" namespace "Fiber":
    cdef cppclass ParallelizationOptions:
        ParallelizationOptions()
        void enableOpenCl(const OpenClOptions&)
        void disableOpenCl()
        cbool isOpenClEnabled() const
        void setMaxThreadCount(int)
        int maxThreadCount() const


cdef class Options:
    cdef:
        AcaOptions aca_options
        ParallelizationOptions parallelization
