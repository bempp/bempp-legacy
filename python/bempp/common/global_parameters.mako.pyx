from cython.operator cimport dereference as deref
from bempp.utils.parameter_list cimport c_ParameterList, ParameterList

def global_parameters():
    cdef c_ParameterList params = c_global_parameters()
    cdef ParameterList p = ParameterList()
    deref(p.impl_).assign(params)
    return p
