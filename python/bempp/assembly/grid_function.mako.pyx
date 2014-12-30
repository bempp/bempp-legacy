
<% from data_types import dtypes, compatible_dtypes, ctypes, scalar_cython_type, real_cython_type
%> 

% for pyvalue in dtypes:
from bempp.utils.armadillo cimport armadillo_to_np_${pyvalue}, armadillo_col_to_np_${pyvalue}
% endfor
    
from bempp.space.space cimport Space
from bempp.utils.armadillo cimport Col
from bempp.utils cimport shared_ptr 
from bempp.space.space cimport c_Space, _py_get_space_ptr 
from bempp.utils.parameter_list cimport ParameterList, c_ParameterList 
from bempp.utils.armadillo cimport Mat
from bempp.utils cimport catch_exception
from bempp.utils cimport complex_float,complex_double
from bempp.utils.enum_types cimport construction_mode
from cython.operator cimport dereference as deref
import numpy as np
cimport numpy as np

np.import_array()

cdef object _fun = None

% for pyvalue,cyvalue in dtypes.items():
cdef void _fun_interface_${pyvalue}(const Col[${real_cython_type(cyvalue)}]& x, 
        const Col[${real_cython_type(cyvalue)}]& normal,
        Col[${cyvalue}]& result):

    cdef np.npy_intp shape_x[1]
    cdef np.npy_intp shape_normal[1]
    cdef np.npy_intp shape_res[1]
    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=1,mode='fortran'] py_x
    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=1,mode='fortran'] py_normal
    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=1,mode='fortran'] py_res
    cdef Col[${real_cython_type(cyvalue)}] my_x = Col[${real_cython_type(cyvalue)}](x)
    cdef Col[${real_cython_type(cyvalue)}] my_normal = Col[${real_cython_type(cyvalue)}](normal)

    shape_x[0] = <np.npy_intp> x.n_rows
    shape_normal[0] = <np.npy_intp> normal.n_rows
    shape_res[0] = <np.npy_intp> result.n_rows

    py_x = np.PyArray_SimpleNewFromData(1,shape_x,np.dtype("${pyvalue}").num,my_x.memptr())
    py_normal = np.PyArray_SimpleNewFromData(1,shape_normal,np.dtype("${pyvalue}").num,my_normal.memptr())
    py_res = np.PyArray_SimpleNewFromData(1,shape_res,np.dtype("${pyvalue}").num,result.memptr())

    _fun(py_x,py_normal,py_res) 
% endfor



cdef class GridFunction:
    
    def __cinit__(self,ParameterList parameter_list,Space space,**kwargs):
        pass

    def __init__(self,ParameterList parameter_list,Space space,**kwargs):

% for pyvalue,cyvalue in dtypes.items():
        cdef Col[${cyvalue}]* arma_data_${pyvalue}
        cdef ${scalar_cython_type(cyvalue)} [:] data_view_${pyvalue}
% endfor

        global _fun

        self._space = space

        if 'basis_type' in kwargs:
            self._basis_type = np.dtype(kwargs['basis_type'])
        else:
            self._basis_type = np.dtype('float64')

        if 'result_type' in kwargs:
            self._result_type = np.dtype(kwargs['result_type'])
        else:
            self._result_type = np.dtype('float64')

        if 'fun' in kwargs:

            _fun = kwargs['fun']

            if 'dual_space' not in kwargs:
                raise ValueError('Need to specify dual space')

            approx_mode = bytes('approximate',"UTF-8")
            if 'approximation_mode' in kwargs:
                approx_mode = bytes(kwargs['approximation_mode'],"UTF-8")

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                self._impl_${pybasis}_${pyresult}.reset(
                        new c_GridFunction[${cybasis},${cyresult}](deref(parameter_list.impl_),
                        _py_get_space_ptr[${cybasis}](self._space.impl_),
                        _py_get_space_ptr[${cybasis}]((<Space>kwargs['dual_space']).impl_),
                        deref(_py_surface_normal_dependent_function_${pyresult}(_fun_interface_${pyresult},3,
                            self._space.codomainDimension)),
                        construction_mode(approx_mode)))
%         endif
%     endfor
% endfor

        elif 'projections' in kwargs:

            if 'dual_space' not in kwargs:
                raise ValueError('Need to specify dual space')

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                num_entries = kwargs['projections'].shape[0]
                data_view_${pyresult} = kwargs['projections']
                arma_data_${pyresult} = new Col[${cyresult}](<${cyresult}*>&data_view_${pyresult}[0],num_entries,True,False)

                self._impl_${pybasis}_${pyresult}.reset(
                        new c_GridFunction[${cybasis},${cyresult}](deref(parameter_list.impl_),
                        _py_get_space_ptr[${cybasis}](self._space.impl_),
                        _py_get_space_ptr[${cybasis}]((<Space>kwargs['dual_space']).impl_),
                        deref(arma_data_${pyresult})))
                del arma_data_${pyresult}
%         endif
%     endfor
% endfor

        elif 'coefficients' in kwargs:

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                num_entries = kwargs['coefficients'].shape[0]
                data_view_${pyresult} = kwargs['coefficients']
                arma_data_${pyresult} = new Col[${cyresult}](<${cyresult}*>&data_view_${pyresult}[0],num_entries,True,False)

                self._impl_${pybasis}_${pyresult}.reset(
                        new c_GridFunction[${cybasis},${cyresult}](deref(parameter_list.impl_),
                        _py_get_space_ptr[${cybasis}](self._space.impl_),
                        deref(arma_data_${pyresult})))
                del arma_data_${pyresult}
%         endif
%     endfor
% endfor
        else:
            raise ValueError("Need to specify at least one of 'coefficients', 'projections' or 'fun'")





% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef np.ndarray _get_coefficients_${pybasis}_${pyresult}(self):
        cdef const Col[${cyresult}]* arma_coeffs = &deref(self._impl_${pybasis}_${pyresult}).coefficients()
        return armadillo_col_to_np_${pyresult}(deref(arma_coeffs))
%          endif
%      endfor
%  endfor



% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef void _set_coefficients_${pybasis}_${pyresult}(self, 
            np.ndarray[${scalar_cython_type(cyresult)},ndim=1,mode='fortran'] coeffs):
        cdef Col[${cyresult}]* arma_coeffs = new Col[${cyresult}](
                <${cyresult}*>&coeffs[0],coeffs.shape[0],False,True)
        deref(self._impl_${pybasis}_${pyresult}).setCoefficients((deref(arma_coeffs)))
        del arma_coeffs
%          endif
%      endfor
%  endfor


% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef np.ndarray _projections_${pybasis}_${pyresult}(self, Space dual_space):
        cdef Col[${cyresult}] arma_coeffs = deref(self._impl_${pybasis}_${pyresult}).projections(
                _py_get_space_ptr[${cybasis}](dual_space.impl_))
        return armadillo_col_to_np_${pyresult}(arma_coeffs)
%          endif
%      endfor
%  endfor


    def projections(self, Space dual_space):

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                return self._projections_${pybasis}_${pyresult}(dual_space)
%          endif
%      endfor
%  endfor

    property coefficients:
        
        def __get__(self):
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                return self._get_coefficients_${pybasis}_${pyresult}()
%          endif
%      endfor
%  endfor

        def __set__(self,coeffs):
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                self._set_coefficients_${pybasis}_${pyresult}(coeffs)
%          endif
%      endfor
%  endfor

    property l2_norm:

        def __get__(self):
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                return deref(self._impl_${pybasis}_${pyresult}).L2Norm_${pybasis}_${pyresult}()
%          endif
%      endfor
%  endfor

    property grid:
        def __get__(self):
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                return self._space.grid
%          endif
%      endfor
%  endfor


    property space:
        def __get__(self):
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                return self._space
%          endif
%      endfor
%  endfor
