from bempp.utils cimport shared_ptr
from libcpp.string cimport string 
from libcpp cimport bool as cbool
from libcpp.vector cimport vector


cdef extern from "Teuchos_ParameterList.hpp" namespace "Teuchos::ParameterList":

    cdef cppclass c_ParameterList "Teuchos::ParameterList":
        c_ParameterList()
        void dump "print" () const
        cbool isParameter(const string& name) const
        cbool isSublist(const string& name) const

        cbool isInt "isType<int>" (const string& name) const
        cbool isDouble "isType<double>" (const string& name) const
        cbool isString "isType<std::string>" (const string& name) const

        int& getInt "get<int>" (const string& name) 
        double& getDouble "get<double>" (const string& name) 
        string& getString "get<std::string>" (const string& name) 
        c_ParameterList& sublist(const string& name) 

        c_ParameterList& setName(const string& name)
        c_ParameterList& setInt "set<int>" (const string& name, const int& value)
        c_ParameterList& setDouble "set<double>" (const string& name, const double& value)
        c_ParameterList& setString "set<std::string>" (const string& name, const string& value)
        c_ParameterList& setList "set<Teuchos::ParameterList>" (const string& name, const c_ParameterList& value)

        cbool remove(const string& name)
        c_ParameterList& setParameters(const c_ParameterList& source)
        int numParams() const

cdef extern from "bempp/utils/utils.hpp" namespace "Bempp":
    int parameter_list_length(const c_ParameterList& parameters)
    vector[string]  parameter_names(const c_ParameterList& parameters)
    string print_parameters(const c_ParameterList& parameters)

cdef class ParameterList:
    cdef c_ParameterList* impl_
    cdef cbool isView_
    cdef int c_len(self)
