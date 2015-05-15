<% from data_types import dtypes, compatible_dtypes, ctypes, scalar_cython_type , ctypes, real_cython_type
%> 
    
from bempp.utils cimport Vector
from bempp.utils cimport shared_ptr 
from bempp.space.space cimport c_Space, _py_get_space_ptr,Space
from bempp.utils.parameter_list cimport ParameterList, c_ParameterList 
from bempp.utils cimport catch_exception
from bempp.utils cimport complex_float,complex_double
from bempp.utils.enum_types cimport ConstructionMode
from bempp.grid.codim_template cimport codim_zero
from bempp.grid.entity cimport Entity0, c_Entity
from bempp.utils.eigen cimport Matrix

cimport numpy as np
import numpy as np


cdef extern from "bempp/fiber/function.hpp" namespace "Fiber":
    cdef cppclass c_Function "Fiber::Function"[RESULT]:
        pass

cdef extern from "bempp/assembly/grid_function.hpp" namespace "Bempp": 
    cdef cppclass c_GridFunction "Bempp::GridFunction"[ BASIS, RESULT ]: 
        c_GridFunction(const c_ParameterList &parameterList,
                       const shared_ptr[c_Space[BASIS]] & space,
                       const Vector[RESULT] & coefficients) except+catch_exception

        c_GridFunction(const c_ParameterList &parameterList,
                       const shared_ptr[c_Space[BASIS]] & space,
                       const shared_ptr[c_Space[BASIS]] & dual_space,
                       const Vector[RESULT] & projections) except+catch_exception

        c_GridFunction(const c_ParameterList &parameterList,
                        const shared_ptr[c_Space[BASIS]]& space,
                        const shared_ptr[c_Space[BASIS]]& dualSpace,
                        const c_Function[RESULT]& function,
                        ConstructionMode constructionMde) except+catch_exception

 
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:

        ${real_cython_type(cyresult)} L2Norm_${pybasis}_${pyresult} "Bempp::GridFunction<${cybasis},${cyresult}>::L2Norm"()

%         endif
%     endfor
% endfor

        const Vector[RESULT]& coefficients() except+catch_exception
        void setCoefficients(const Vector[RESULT]& coeffs) except+catch_exception
        Vector[RESULT] projections(const shared_ptr[const c_Space[BASIS]] &dualSpace)


        void evaluate(const c_Entity[codim_zero]& element,
                const Matrix[BASIS]& local,
                Matrix[RESULT]& values)

    
cdef extern from "bempp/assembly/py_functors.hpp" namespace "Bempp":
% for pyvalue,cyvalue in dtypes.items():
    cdef shared_ptr[c_Function[${cyvalue}]] _py_surface_normal_dependent_function_${pyvalue} "Bempp::_py_surface_normal_dependent_function<${ctypes(cyvalue)}>"(
            void (*callable)(object,object,int, object, object),object,
            int argumentDimension, int resultDimension) except+catch_exception
% endfor

cdef class GridFunction:
    cdef object _basis_type,
    cdef object _result_type,
    cdef Space _space,
    cdef ParameterList _parameter_list

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef shared_ptr[c_GridFunction[${cybasis},${cyresult}]] _impl_${pybasis}_${pyresult}

    cdef np.ndarray _get_coefficients_${pybasis}_${pyresult}(self)
    cdef void _set_coefficients_${pybasis}_${pyresult}(self,np.ndarray[${scalar_cython_type(cyresult)},ndim=1,mode='fortran'] coeffs)
    cdef np.ndarray _projections_${pybasis}_${pyresult}(self,Space dual_space)
    cdef np.ndarray _evaluate_complex(self, Entity0 element, np.ndarray local_coordinates)
    cdef np.ndarray _evaluate_real(self, Entity0 element, np.ndarray local_coordinates)
%         endif
%     endfor
% endfor    



