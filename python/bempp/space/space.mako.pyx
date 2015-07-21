""" Wrappers for all types of C++ spaces """
<%
from data_types import dtypes, real_cython_type
from space import spaces
ifloop = lambda x: 'if' if loop.index == 0 else 'elif'
%>
from bempp.grid.grid cimport Grid, c_Grid
from cython.operator cimport dereference as deref
from libcpp cimport bool as cbool
from bempp.utils.eigen cimport eigen_matrix_to_np_float32,eigen_matrix_to_np_float64
from bempp.utils.shared_ptr cimport const_pointer_cast

cdef class Space:
    """ Space of functions defined on a grid

        The exact space depends on the input.
    """
    def __init__(self, Grid grid not None, unsigned int order, object type_string):
        super(Space, self).__init__()
        self._order = order
        self._type_string = type_string
        self._grid = grid

    def update(self):

        from bempp import function_space
        cdef Space tmp_space = function_space(self._grid,self._type_string,self._order)
        self.impl_.assign(tmp_space.impl_)

    property dtype:
        """ Type of the basis functions in this space. """
        def __get__(self):
            from numpy import dtype
            return dtype(self.impl_.dtype());

    property codomain_dimension:
        """ Number of components of values of functions in this space (e.g. 1 for scalar functions). """
        def __get__(self):
            return self.impl_.codomainDimension()

    property grid:
        """ The underlyign grid for the space. """
        def __get__(self):
            cdef Grid result = Grid.__new__(Grid)
            result.impl_ = const_pointer_cast[c_Grid](self.impl_.grid())
            return result

    property domain_dimension:
        """ Dimension of the domain on which the space is defined (default 2). """

        def __get__(self):
            return self.impl_.domainDimension()

    property global_dof_count:
        """ Number of global degrees of freedom. """

        def __get__(self):
            return self.impl_.globalDofCount()

    property flat_local_dof_count:
        """ Total number of local degrees of freedom. """

        def __get__(self):
            return self.impl_.flatLocalDofCount()

    cpdef cbool is_compatible(self,Space other):
        """ s.is_compatible(other_space)

            Test if both spaces have the same global degrees of freedom.

        """


        if not (self.dtype==other.dtype): return False
        return self.impl_.isCompatible(other.impl_)

    property order:
        """ The order of the basis functions. """

        def __get__(self):
            return self._order


    def get_global_dofs(self,Entity0 element):

        cdef vector[int] global_dofs_vec
        cdef vector[double] local_dof_weights_vec
        if self.dtype=='float64':
            _py_space_get_global_dofs_float64(self.impl_,deref(element.impl_),global_dofs_vec,local_dof_weights_vec)
        else:
            raise ValueError("Unsupported dtype.")
        global_dofs = []
        local_dof_weights = []
        #for i in range(global_dofs_vec.size()):
        #    global_dofs.append(global_dofs_vec[i])
        #for i in range(local_dof_weights_vec):
        #    local_dof_weights.append(local_dof_weights_vec[i])
        #return (global_dofs,local_dof_weights)
        return (global_dofs_vec,local_dof_weights_vec) 
        
           

                

% for class_name, description in spaces.items():
cdef class ${class_name}(Space):
    """ ${description['doc']}

        Attributes
        ----------

        grid : Grid
            Grid over which to discretize the space.

        dtype : numpy.dtype
            Type of the basis functions in this space.

        codomain_dimension : int
            Number of components of values of functions in this space.

        domain_dimension : int
            Dimension of the domain on which the space is defined.

        global_dof_count : int
            Number of global degrees of freedom.

        flat_local_dof_count : int
            Total number of local degrees of freedom.

        global_dof_interpolation_points : np.ndarray 
            (3xN) matrix of global interpolation points for the space,
            where each column is the coordinate of an interpolation point.

        global_dof_interpolation_points : np.ndarray
            (3xN) matrix of normal directions associated with the interpolation points.

        order : int
            Order of the polynomial degre of the space. 

        Notes
        -----
        A space instance should always be created using the function 'bempp.function_space'.
    """
    def __init__(self, Grid grid not None, order, object type_string):

        from numpy import dtype as np_dtype
        super(${class_name}, self).__init__(grid,order,type_string)

        dtype = 'float64'
        if dtype not in ${list(dtypes.keys())}:
                raise TypeError("Unexpected basis type")
%    for pytype, cytype in dtypes.items():
        if dtype == "${pytype}":
%       if description['implementation'] == 'grid_only':
            self.impl_.set( shared_ptr[c_Space[${cytype}]](
                <c_Space[${cytype}]*>
                new ${'c_' + class_name}[${cytype}](grid.impl_)
            ))
%       elif description['implementation'] == 'polynomial':
            self.impl_.set( shared_ptr[c_Space[${cytype}]](
                <c_Space[${cytype}]*> new ${'c_' + class_name}[${cytype}](
                    grid.impl_, <int> self.order
                )
            ))
%       endif
%    endfor


    property global_dof_interpolation_points:
        """ (3xN) matrix of global interpolation points for the space, where each column is the
            coordinate of an interpolation points. """

        def __get__(self):

% for pybasis,cybasis in dtypes.items():
            if self.dtype=="${pybasis}":
%     if pybasis in ['float32','complex64']:
                    return eigen_matrix_to_np_float32(_py_space_get_global_dof_interp_points_${pybasis}(self.impl_))
%     else:
                    return eigen_matrix_to_np_float64(_py_space_get_global_dof_interp_points_${pybasis}(self.impl_))
%     endif
% endfor
            raise("Unknown dtype for space")


    property global_dof_normals:
        """ (3xN) matrix of normal directions associated with the interpolation points. """

        def __get__(self):

% for pybasis,cybasis in dtypes.items():
            if self.dtype=="${pybasis}":
%     if pybasis in ['float32','complex64']:
                    return eigen_matrix_to_np_float32(_py_space_get_global_dof_normals_${pybasis}(self.impl_))
%     else:
                    return eigen_matrix_to_np_float64(_py_space_get_global_dof_normals_${pybasis}(self.impl_))
%     endif
% endfor
            raise("Unknown dtype for space")


% endfor



cdef class RaviartThomas0VectorSpace(Space):
    """ Raviart Thomas Basis functions of order 0

        Attributes
        ----------

        grid : Grid
            Grid over which to discretize the space.

        dtype : numpy.dtype
            Type of the basis functions in this space.

        codomain_dimension : int
            Number of components of values of functions in this space.

        domain_dimension : int
            Dimension of the domain on which the space is defined.

        global_dof_count : int
            Number of global degrees of freedom.

        flat_local_dof_count : int
            Total number of local degrees of freedom.

        global_dof_interpolation_points : np.ndarray 
            (3xN) matrix of global interpolation points for the space,
            where each column is the coordinate of an interpolation point.

        global_dof_interpolation_points : np.ndarray
            (3xN) matrix of normal directions associated with the interpolation points.

        order : int
            Order of the polynomial degre of the space. 

        Notes
        -----
        A space instance should always be created using the function 'bempp.function_space'.
    """

    def __init__(self, Grid grid not None, order, object type_string):

        super(RaviartThomas0VectorSpace,self).__init__(grid,order, type_string)
        self.impl_.set( shared_ptr[c_Space[double]](
            <c_Space[double]*>new c_RaviartThomas0VectorSpace(grid.impl_)))

