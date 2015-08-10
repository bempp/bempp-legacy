#cython: embedsignature=True
<% from data_types import dtypes, compatible_dtypes, ctypes, scalar_cython_type, real_cython_type
%> 

% for pyvalue in dtypes:
from bempp.utils.eigen cimport eigen_vector_to_np_${pyvalue}
from bempp.utils.eigen cimport np_to_eigen_vector_${pyvalue}
from bempp.utils.eigen cimport np_to_eigen_matrix_${pyvalue}
from bempp.utils.eigen cimport eigen_matrix_to_np_${pyvalue}
% endfor

from bempp.utils.eigen cimport Matrix
from bempp.space.space cimport Space
from bempp.utils cimport Matrix,Vector
from bempp.utils cimport shared_ptr 
from bempp.space.space cimport c_Space, _py_get_space_ptr 
from bempp.utils.parameter_list cimport ParameterList, c_ParameterList 
from bempp.utils cimport catch_exception
from bempp.utils cimport complex_float,complex_double
from bempp.utils.enum_types cimport construction_mode
from bempp.grid.entity cimport Entity0
from bempp.common import global_parameters
from cython.operator cimport dereference as deref
import numpy as np
cimport numpy as np


cdef void _fun_interface(object x, object normal, int domain_index, object res, object call_fun):

    call_fun(x,normal,domain_index,res) 


cdef class GridFunction:
    """

    This class represents functions defined on a grid. It can be initialized
    in three different ways.

    1. By providing a Python callable. Any Python callable of the following form
       is valid.::

            callable(x,n,domain_index,result)

       Here, x, n, and result are all numpy arrays. x contains the current evaluation
       point, n the associated outward normal direction and result is a numpy array
       that will store the result of the Python callable. The variable domain_index
       stores the index of the subdomain on which x lies (default 0). This makes it
       possible to define different functions for different subdomains.

       The following example defines input data that is the inner product of the
       coordinate x with the normal direction n.::

            fun(x,n,domain_index,result):
                result[0] =  np.dot(x,n)

       If the input function returns complex data the keyword argument 
       'complex_data=True' needs to be specified in the contructor of the GridFunction.    
    2. By providing a vector of coefficients at the nodes. This is preferable if
       the coefficients of the data are coming from an external code.

    3. By providing a vector of projection data and a corresponding dual space.


    Parameters
    ----------
    space : bempp.Space
        The space over which the GridFunction is defined.
    complex_data : bool
        Specify whether an input function returns complex numbers. 
        (optional, default False).
    dual_space : bempp.Space
        A representation of the dual space (optional).
    fun : callable
        A Python function from which the GridFunction is constructed
        (optional).
    coefficients : np.ndarray
        A 1-dimensional array with the coefficients of the GridFunction
        at the interpolatoin points of the space (optional).
    projections : np.ndarray
        A 1-dimensional array with the projections of the GridFunction
        onto a dual space (optional).
    parameter_list : bempp.ParameterList
        A ParameterList object used for the assembly of
        the GridFunction (optional).

    Attributes
    ----------
    coefficients : np.ndarray
        Return or set the vector of coefficients.
    l2_norm : double
        Return the L2 norm of the GridFunction.
    space : bemp.Space
        Return the space over which the GridFunction is defined.
    grid : bempp.Grid
        Return the underlying grid.
    parameter_list : bempp.ParameterList
        Return the set of parameters.
    basis_type : np.dtype
        Return the basis function type.
    result_type : np.dtype
        Return the result type.


    Notes
    -----
    * Only one of projections, coefficients, or fun is allowed as parameter.
    * To export a GridFunction to a file see the module bempp.file_interfaces.

    Examples
    --------
    To create a GridFunction from a real Python callable my_fun use

    >>> grid_function = GridFunction(space, fun=my_fun)

    To create a GridFunction from a complex Python callable my_fun use

    >>> grid_function = GridFunction(space, fun=my_fun,
    ...    complex_data=True)

    To create a GridFunction from a vector of coefficients coeffs use

    >>> grid_function = GridFunction(space,coefficients=coeffs)

    To create a GridFunction from a vector of projections proj use
    
    >>> grid_function = GridFunction(space,dual_space=dual_space, projections=proj)


    """

    
    def __cinit__(self,Space space,**kwargs):
        pass

    def __init__(self,Space space,**kwargs):

% for pyvalue,cyvalue in dtypes.items():
        cdef Vector[${cyvalue}]* eigen_data_${pyvalue}
        cdef ${scalar_cython_type(cyvalue)} [::1] data_view_${pyvalue}
% endfor

        if 'parameter_list' in kwargs:
            self._parameter_list = kwargs['parameter_list']
        else:
            import bempp
            self._parameter_list = bempp.global_parameters

        global _fun

        self._space = space

        self._basis_type = space.dtype
        self._result_type = None


        if 'fun' in kwargs:

            if 'complex_data' in kwargs:
                self._result_type = np.dtype('complex128')
            else:
                self._result_type = np.dtype('float64')

            if 'dual_space' not in kwargs:
                dual_space = space
            else:
                dual_space = kwargs['dual_space']

            approx_mode = 'approximate'.encode("UTF-8")
            if 'approximation_mode' in kwargs:
                approx_mode = kwargs['approximation_mode'].encode("UTF-8")

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                self._impl_${pybasis}_${pyresult}.reset(
                        new c_GridFunction[${cybasis},${cyresult}](deref((<ParameterList>self.parameter_list).impl_),
                        _py_get_space_ptr[${cybasis}](self._space.impl_),
                        _py_get_space_ptr[${cybasis}]((<Space>dual_space).impl_),
                        deref(_py_surface_normal_dependent_function_${pyresult}(_fun_interface,kwargs['fun'],3,
                            self._space.codomain_dimension)),
                        construction_mode(approx_mode)))
%         endif
%     endfor
% endfor

        elif 'projections' in kwargs:

            if np.iscomplexobj(kwargs['projections']):
                self._result_type = np.dtype('complex128')
            else:
                self._result_type = np.dtype('float64')

            if 'dual_space' not in kwargs:
                raise ValueError('Need to specify dual space')

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                num_entries = kwargs['projections'].shape[0]
                eigen_data_${pyresult} = new Vector[${cyresult}](np_to_eigen_vector_${pyresult}(np.require(kwargs['projections'],"${pyresult}","F")))

                self._impl_${pybasis}_${pyresult}.reset(
                        new c_GridFunction[${cybasis},${cyresult}](deref((<ParameterList>self.parameter_list).impl_),
                        _py_get_space_ptr[${cybasis}](self._space.impl_),
                        _py_get_space_ptr[${cybasis}]((<Space>kwargs['dual_space']).impl_),
                        deref(eigen_data_${pyresult})))
                del eigen_data_${pyresult}
%         endif
%     endfor
% endfor

        elif 'coefficients' in kwargs:

            if np.iscomplexobj(kwargs['coefficients']):
                self._result_type = np.dtype('complex128')
            else:
                self._result_type = np.dtype('float64')

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                num_entries = kwargs['coefficients'].shape[0]
                eigen_data_${pyresult} = new Vector[${cyresult}](np_to_eigen_vector_${pyresult}(np.require(kwargs['coefficients'],"${pyresult}","F")))


                self._impl_${pybasis}_${pyresult}.reset(
                        new c_GridFunction[${cybasis},${cyresult}](deref((<ParameterList>self.parameter_list).impl_),
                        _py_get_space_ptr[${cybasis}](self._space.impl_),
                        deref(eigen_data_${pyresult})))
                del eigen_data_${pyresult}
%         endif
%     endfor
% endfor
        else:
            raise ValueError("Need to specify at least one of 'coefficients', 'projections' or 'fun'")

    def __dealloc__(self):


% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
        self._impl_${pybasis}_${pyresult}.reset()
%         endif
%     endfor
% endfor    



% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef np.ndarray _get_coefficients_${pybasis}_${pyresult}(self):
        cdef const Vector[${cyresult}]* eigen_coeffs = &deref(self._impl_${pybasis}_${pyresult}).coefficients()
        return eigen_vector_to_np_${pyresult}(deref(eigen_coeffs))
%          endif
%      endfor
%  endfor



% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef void _set_coefficients_${pybasis}_${pyresult}(self, 
            np.ndarray coeffs):
        deref(self._impl_${pybasis}_${pyresult}).setCoefficients(np_to_eigen_vector_${pyresult}(np.require(coeffs,"${pyresult}","F")))
%          endif
%      endfor
%  endfor



% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
    cdef np.ndarray _projections_${pybasis}_${pyresult}(self, Space dual_space):
        cdef Vector[${cyresult}] eigen_coeffs = deref(self._impl_${pybasis}_${pyresult}).projections(
                _py_get_space_ptr[${cybasis}](dual_space.impl_))
        return eigen_vector_to_np_${pyresult}(eigen_coeffs)
%          endif
%      endfor
%  endfor

    cdef np.ndarray _evaluate_complex(self, Entity0 element, np.ndarray local_coordinates):
        cdef Matrix[complex_double] values
        deref(self._impl_float64_complex128).evaluate(
                deref(element.impl_),np_to_eigen_matrix_float64(local_coordinates),
                values)
        return eigen_matrix_to_np_complex128(values)

    cdef np.ndarray _evaluate_real(self, Entity0 element, np.ndarray local_coordinates):
        cdef Matrix[double] values
        deref(self._impl_float64_float64).evaluate(
                deref(element.impl_),np_to_eigen_matrix_float64(local_coordinates),
                values)
        return eigen_matrix_to_np_float64(values)


    def evaluate(self, element, local_coordinates):

        if self._basis_type!='float64':
            raise ValueError("Only basis_type 'float64' is supported.")
        if self._result_type=='complex128':
            return self._evaluate_complex(element, local_coordinates)
        elif self._result_type=='float64':
            return self._evaluate_real(element, local_coordinates)
        else:
            raise ValueError("Result type {0} not supported.".format(str(self._result_type)))


    def projections(self, Space dual_space):
        """

        Compute the vector of projections onto the 
        given dual space.

        Parameters
        ----------
        dual_space : bempp.Space
            A representation of the dual space.

        Returns
        -------
        out : np.ndarray
            A vector of projections onto the dual space.

        """

% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
        if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
            return self._projections_${pybasis}_${pyresult}(dual_space)
%          endif
%      endfor
%  endfor


    def plot(self):
        """

        Plot a grid function using Gmsh.

        """

        from bempp.external.viewers import visualize_with_gmsh
        visualize_with_gmsh(self)


    def __add__(self,GridFunction other):

        if not self.space.is_compatible(other.space):
            raise ValueError("Spaces do not match")

        return GridFunction(self.space,
                coefficients=self.coefficients+other.coefficients,
                parameter_list=self.parameter_list)


    def __mul__(self,object alpha):

        if not isinstance(self,GridFunction):
            return alpha*self

        if np.isscalar(alpha):
            return GridFunction(self.space,
                    coefficients=alpha*self.coefficients,
                    parameter_list=self.parameter_list)
        else:
            raise NotImplementedError("Cannot multiply Gridfunction with object of type "+str(type(alpha)))

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self,GridFunction other):
        return self.__add__(-other)

    property coefficients:
        """ Return or set the vector of coefficients. """

        
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
        """ Return the L2 norm of the GridFunction. """

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
        """ Return the underlying grid. """

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
        """ Return the space over which the GridFunction is defined. """

        def __get__(self):
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:
            if (self._basis_type=="${pybasis}") and (self._result_type=="${pyresult}"):
                return self._space
%          endif
%      endfor
%  endfor

    property parameter_list:
        """ Return the set of parameters. """

        def __get__(self):
            return self._parameter_list

    property basis_type:
        """ Return the basis function type. """

        def __get__(self):
            return self._basis_type

    property result_type:
        """ Return the result type. """

        def __get__(self):
            return self._result_type
