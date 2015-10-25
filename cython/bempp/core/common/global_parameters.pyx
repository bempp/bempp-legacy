from bempp.core.utils cimport c_ParameterList, ParameterList
from cython.operator cimport dereference as deref

cdef extern from "bempp/common/global_parameters.hpp" namespace "Bempp":
    c_ParameterList c_global_parameters "Bempp::GlobalParameters::parameterList" ()

def global_parameters():
    cdef c_ParameterList params = c_global_parameters()
    cdef ParameterList p = ParameterList()
    deref(p.impl_).assign(params)
    return p
