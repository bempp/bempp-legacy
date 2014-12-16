from bempp.utils.parameter_list cimport c_ParameterList

cdef extern from "common/global_parameters.hpp" namespace "Bempp":
    c_ParameterList c_global_parameters "Bempp::GlobalParameters::parameterList" ()




