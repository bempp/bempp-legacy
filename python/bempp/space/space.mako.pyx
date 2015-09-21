""" Wrappers for all types of C++ spaces """
<%
from data_types import dtypes, real_cython_type
from space import spaces
ifloop = lambda x: 'if' if loop.index == 0 else 'elif'
%>
from bempp.grid.grid cimport Grid
from cython.operator cimport dereference as deref
from libcpp cimport bool as cbool
from bempp.utils.shared_ptr cimport const_pointer_cast
from bempp.utils.eigen cimport eigen_matrix_to_np_float32,eigen_matrix_to_np_float64



cdef class Space:
    """ Space of functions defined on a grid

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

        global_dof_interpolation_normals : np.ndarray
            (3xN) matrix of normal directions associated with the interpolation points.

        order : int
            Order of the polynomial degre of the space. 

        Notes
        -----
        A space instance should always be created using the function 'bempp.function_space'.

    """
    def __init__(self, unsigned int order, str kind, domains=None,
                 cbool closed=True, gridname=None):
        super(Space, self).__init__()
        self._order = order
        self._kind = kind
        self._domains = domains
        self._closed = closed
        self._gridname = gridname

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

    property kind:
        """ The type of space. """

        def __get__(self):
            return self._kind

    property closed:
        """ Specifies whether the space is defined on a closed
        or open subspace. """

        def __get__(self):
            return self._closed

    property domains:
        """ List of integers specifying a list of physical entities
        of subdomains that should be included in the space. None
        is all domains. """

        def __get__(self):
            return self._domains

    property gridname:
        """ Name of the variable that holds the grid used to create the space.
        This is only used to serialize the space with a reference to the grid"""

        def __get__(self):
            return self._gridname

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


    def __getstate__(self):
        state = dict()
        state['kind'] =  self.kind
        state['order'] = self.order
        state['domains'] = self.domains
        state['closed'] = self.closed

        return state

    def __setstate__(self, state):
        raise NotImplementedError('Cannot serialize Space alone. Please '
                                  'use a bempp.utils.ParallelInterface to '
                                  'serialize a spaces with matching grid')


def function_space(Grid grid, kind, order, domains=None, cbool closed=True,
                   gridname='none'):
    """ 

    Return a space defined over a given grid.

    Parameters
    ----------
    grid : bempp.Grid
        The grid object over which the space is defined.

    kind : string
        The type of space. Currently, the following types
        are supported:
        "P" : Continuous and piecewise polynomial functions.
        "DP" : Discontinuous and elementwise polynomial functions.

    order : int
        The order of the space, e.g. 0 for piecewise const, 1 for
        piecewise linear functions.

    domains : list
        List of integers specifying a list of physical entities
        of subdomains that should be included in the space. None
        is all domains.

    closed : bool
        Specifies whether the space is defined on a closed
        or open subspace.

    gridname: str
        name of the variable used to hold the grid. Used only for IPython
        serialization. This is ugly and should be changed.

    Notes
    -----
    The most frequent used types are the space of piecewise constant
    functions (kind="DP", order=0) and the space of continuous,
    piecewise linear functions (kind="P", order=1).

    This is a factory function that initializes a space object. To 
    see a detailed help for space objects see the documentation
    of the instantiated object.

    Examples
    --------
    To initialize a space of piecewise constant functions use

    >>> space = function_space(grid,"DP",0)

    To initialize a space of continuous, piecewise linear functions, use

    >>> space = function_space(grid,"P",1)

    """

    cdef Space s = Space(order, kind, domains, closed, gridname)
    __function_space(s, grid, kind, order, domains, closed)
    return s


def __function_space(Space s, Grid grid, kind, order, domains=None, cbool closed=True):
    """
        helper function for function_space and unpickleing of spaces.
    """
    if kind=="P":
        if not (order>=1 and order <=10):
            raise ValueError("Order must be between 1 and 10")
        if (order==1):
            if domains is None:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpace[double](grid.impl_)))
            else:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpace[double](grid.impl_, domains, closed)))
        else:
            if domains is None:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewisePolynomialContinuousScalarSpace[double](grid.impl_, order)))
            else:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewisePolynomialContinuousScalarSpace[double](grid.impl_, order, domains, closed)))
    elif kind=="DP":
        if not (order>=0 and order <=10):
            raise ValueError("Order must be between 0 and 10")
        if (order==0):
            if domains is None:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpace[double](grid.impl_)))
            else:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpace[double](grid.impl_, domains, closed)))
        else:
            if domains is None:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewisePolynomialDiscontinuousScalarSpace[double](grid.impl_, order)))
            else:
                s.impl_.set(shared_ptr[c_Space[double]](adaptivePiecewisePolynomialDiscontinuousScalarSpace[double](grid.impl_, order, domains, closed)))
    elif kind=="RT":
        if order!=0:
            raise ValueError("Only 0 order Raviart-Thomas spaces are implemented.")
        if domains is None:
            s.impl_.set(shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpace[double](grid.impl_)))
        else:
            s.impl_.set(shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpace[double](grid.impl_, domains, closed)))
    else:
        raise ValueError("Unknown kind")

    return s
