from cython.operator cimport dereference as deref
from cython.operator cimport address
from libcpp.string cimport string 
from libcpp cimport bool as cbool
from .byte_conversion import convert_to_bytes
from libcpp.vector cimport vector

cdef class ParameterList:

    def __cinit__(self):
        self.impl_ = new c_ParameterList()
        self.isView_ = False

    def __dealloc__(self):
        if not self.isView_: del self.impl_

    def _is_parameter(self, name):
        s = convert_to_bytes(name)
        return deref(self.impl_).isParameter(s)

    def _type(self, name):
        s = convert_to_bytes(name)

        if not self._is_parameter(s):
            raise KeyError(name)

        if deref(self.impl_).isSublist(s):
            return 'l'

        if deref(self.impl_).isInt(s):
            return 'i'

        if deref(self.impl_).isDouble(s):
            return 'd'

        if deref(self.impl_).isString(s):
            return 's'

    def _set_parameter(self, name, value):

        s = convert_to_bytes(name)

        if isinstance(value,int):
            deref(self.impl_).setInt(s,value)
        elif isinstance(value,float):
            deref(self.impl_).setDouble(s,value)
        elif isinstance(value,str):
            deref(self.impl_).setString(s,convert_to_bytes(value))
        elif isinstance(value,dict):
            p = ParameterList()
            p.insert_dict(value)
            deref(self.impl_).setList(s,deref(p.impl_))
        else:
            raise ValueError('value must be one of int, float or string.')

    def _get_parameter(self,name):

        cdef string ret_string
        s = convert_to_bytes(name)

        if not self._is_parameter(name):
            raise KeyError(name)

        if self._type(name)=='l':
            return self._sublist(name)
        elif self._type(name)=='i':
            return deref(self.impl_).getInt(s)
        elif self._type(name)=='d':
            return deref(self.impl_).getDouble(s)
        elif self._type(name)=='s':
            ret_string = deref(self.impl_).getString(s)
            return ret_string.decode("UTF-8")
        else:
            raise ValueError(name)

    def _sublist(self,name):
        s = convert_to_bytes(name)
        cdef ParameterList p = ParameterList()
        del p.impl_
        p.impl_ = address(deref(self.impl_).sublist(s))
        p.isView_ = True
        return p

    def set_name(self, name):

        s = convert_to_bytes(name)
        deref(self.impl_).setName(s)

    def _remove(self, name):

        s = convert_to_bytes(name)
        if not self._is_parameter(name):
            raise ValueError(name)
        deref(self.impl_).remove(s)

    def insert_dict(self, d):

        for key in d:
            s = convert_to_bytes(key)
            if isinstance(d[key],dict):
                sublist = self._sublist(s)
                sublist.insert_dict(d[key])
            else:
                self._set_parameter(s,d[key])

    @classmethod
    def from_dict(self,d):
        p = ParameterList()
        p.insert_dict(d)
        return p

    cdef int c_len(self):

        return deref(self.impl_).numParams()
        
    def __len__(self):
        return self.c_len()

    def to_dict(self):

        cdef vector[string] names = parameter_names(deref(self.impl_))
        d = dict()
        d = {name.decode("UTF-8"): self._get_parameter(name).to_dict() if self._type(name) == 'l' else self._get_parameter(name) for name in names}
        return d

    def __getitem__(self,key):

        if not isinstance(key,str):
            raise TypeError(key)

        if not self._is_parameter(key):
            raise KeyError(key)

        return self._get_parameter(key)

    def __setitem__(self,key,value):

        if not isinstance(key,str):
            raise TypeError(key)

        self._set_parameter(key,value)

    def __delitem__(self,key):
        self._remove(key)

    def __iter__(self):

        cdef vector[string] names = parameter_names(deref(self.impl_))

        for s in names:
            yield s.decode("UTF-8")

    def __missing__(self,key):

        return key

    def __contains__(self,key):

        return self._is_parameter(key)

    def __str__(self):
        cdef string s = print_parameters(deref(self.impl_))
        return s.decode("UTF-8")


    




        

        



