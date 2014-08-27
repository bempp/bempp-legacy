from libcpp cimport bool as cbool
<%
from bempp_options import options

def transform(type):
    return {'bool': 'cbool'}.get(type, type)
%>

# Declare enums from various files
cdef extern from "bempp/assembly/aca_options.hpp":
    cdef enum AcaAssemblyMode "Bempp::AcaOptions::AcaAssemblyMode":
        GLOBAL_ASSEMBLY "Bempp::AcaOptions::GLOBAL_ASSEMBLY"
        LOCAL_ASSEMBLY "Bempp::AcaOptions::LOCAL_ASSEMBLY"
        HYBRID_ASSEMBLY "Bempp::AcaOptions::HYBRID_ASSEMBLY"
    cdef enum ReactionToUnsupportedMode \
            "Bempp::AcaOptions::ReactionToUnsupportedMode":
        WARNING "Bempp::AcaOptions::WARNING"
        ERROR "Bempp::AcaOptions::ERROR"
        IGNORE "Bempp::AcaOptions::IGNORE"

# Declare c structure used in conversions
cdef extern from "bempp/assembly/aca_options.hpp" namespace "Bempp":
    cdef cppclass AcaOptions:
        AcaOptions()
% for name, description in options.iteritems():
<%
    if description['c origin'] != 'AcaOptions':
        continue
%>
        ${transform(description['type'])} ${name}
% endfor

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

cdef extern from "bempp/fiber/verbosity_level.hpp" namespace "Fiber":
    cdef enum VerbosityLevel "Bempp::VerbosityLevel::Level":
        VERBOSITY_MEDIUM "Bempp::VerbosityLevel::DEFAULT"
        VERBOSITY_HIGH "Bempp::VerbosityLevel::HIGH"
        VERBOSITY_LOW "Bempp::VerbosityLevel::LOW"

cdef extern from "bempp/assembly/assembly_options.hpp" namespace "Bempp":
    cdef enum BlasQuadrature "Bempp::AssemblyOptions::Value":
        BLAS_QUADRATURE_AUTO "Bempp::AssemblyOptions::AUTO"
        BLAS_QUADRATURE_NO "Bempp::AssemblyOptions::NO"
        BLAS_QUADRATURE_YES "Bempp::AssemblyOptions::YES"
    cdef enum AssemblyMode "Bempp::AssemblyOptions::Mode":
        ASSEMBLY_MODE_DENSE "Bempp::AssemblyOptions::DENSE"
        ASSEMBLY_MODE_ACA "Bempp::AssemblyOptions::ACA"
    cdef cppclass AssemblyOptions:
        AssemblyOptions()
        void switchToDenseMode()
        void switchToAcaMode(const AcaOptions& acaOptions)
        AssemblyMode assemblyMode() const
% for option, description in options.iteritems():
<%
    if description['c origin'] != 'AssemblyOptions':
        continue
%>
        ${description['type']} ${description['getter']}()
        void ${description['setter']}(${description['type']})
% endfor

cdef class Options:
    cdef:
        AcaOptions aca_options
        ParallelizationOptions parallelization
        AssemblyOptions assembly
