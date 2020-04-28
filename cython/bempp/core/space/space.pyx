from cython.operator cimport dereference as deref
from bempp.core.utils cimport Matrix
from bempp.core.utils.shared_ptr cimport reverse_const_pointer_cast
from bempp.core.utils.shared_ptr cimport const_pointer_cast
from bempp.core.utils cimport eigen_matrix_to_np_float64
from bempp.core.utils cimport eigen_matrix_to_np_complex128
from bempp.core.utils cimport np_to_eigen_matrix_float64
from bempp.core.utils cimport np_to_eigen_matrix_complex128
from bempp.core.utils cimport np_to_eigen_vector_float64
from bempp.core.utils cimport np_to_eigen_vector_complex128
from bempp.core.fiber cimport Shapeset
from libcpp.vector cimport vector

cdef class LocalDof:
    cdef int entity_index
    cdef int dof_index

    def __cinit__(self, int entity_index, int dof_index):
        pass

    def __init__(self, int entity_index, int dof_index):
        self.entity_index = entity_index
        self.dof_index = dof_index

    property entity_index:

        def __get__(self):
            return self.entity_index

    property dof_index:

        def __get__(self):
            return self.dof_index

cdef class Space:


    def __cinit__(self):
        pass

    def __init__(self):
        super(Space, self).__init__()

    def __dealloc__(self):
        self.impl_.reset()

    property dtype:
        """ Type of the basis functions in this space. """
        def __get__(self):
            from numpy import dtype
            return dtype('float64');

    property codomain_dimension:
        """ Number of components of values of functions in this space (e.g. 1 for scalar functions). """
        def __get__(self):
            return deref(self.impl_).codomainDimension()

    property grid:
        """ The underlyign grid for the space. """
        def __get__(self):
            cdef Grid result = Grid.__new__(Grid)
            result.impl_ = const_pointer_cast[c_Grid](deref(self.impl_).grid())
            return result

    property domain_dimension:
        """ Dimension of the domain on which the space is defined (default 2). """

        def __get__(self):
            return deref(self.impl_).domainDimension()

    property global_dof_count:
        """ Number of global degrees of freedom. """

        def __get__(self):
            return deref(self.impl_).globalDofCount()

    property flat_local_dof_count:
        """ Total number of local degrees of freedom. """

        def __get__(self):
            return deref(self.impl_).flatLocalDofCount()

    property discontinuous_space:
        """Return the associated discontinuous scalar space."""

        def __get__(self):
            cdef Space space = Space()
            space.impl_.assign(deref(self.impl_).discontinuousSpace(self.impl_))
            return space

    property is_discontinuous:
        """Return true of basis functions are scalar and only extend over a single element."""

        def __get__(self):
            return deref(self.impl_).isDiscontinuous()

    def is_compatible(self,Space other):
        """ Test if both spaces have the same global degrees of freedom. """

        return deref(self.impl_).spaceIsCompatible(deref(other.impl_))

    def __richcmp__(Space self, Space other, int op):
        """Comparison operator."""

        if op == 2:
            return self.is_identical(other)
        elif op == 3:
            return not self.is_identical(other)
        else:
            return NotImplemented

    def is_identical(self, Space other):
        return self.impl_.get() == other.impl_.get()

    def get_global_dofs(self,Entity0 element, dof_weights=False):

        cdef vector[int] global_dofs_vec
        cdef vector[double] local_dof_weights_vec
        deref(self.impl_).getGlobalDofs(deref(element.impl_), global_dofs_vec, local_dof_weights_vec)
        if dof_weights:
            return global_dofs_vec,local_dof_weights_vec
        else:
            return global_dofs_vec

    def shapeset(self, Entity0 element):

        cdef Shapeset shapeset = Shapeset()
        shapeset.impl_ = &deref(self.impl_).shapeset(deref(element.impl_))
        return shapeset


    def evaluate_local_basis(self, Entity0 element, object local_coordinates, object local_coefficients):
        """Evaluate local basis functions on a given element."""

        import numpy as np

        if np.isreal(local_coefficients).all():
            coeffs_real = local_coefficients
            return eigen_matrix_to_np_float64(
                    c_evaluateLocalBasis[double](deref(self.impl_), deref(element.impl_),
                                                np_to_eigen_matrix_float64(local_coordinates),
                                                np_to_eigen_vector_float64(local_coefficients)))
        else:
            coeffs_complex = local_coefficients
            return eigen_matrix_to_np_complex128(
                    c_evaluateLocalBasis[complex_double](deref(self.impl_), deref(element.impl_),
                                                        np_to_eigen_matrix_float64(local_coordinates),
                                                        np_to_eigen_vector_complex128(local_coefficients)))


    def global_to_local_dofs(self, global_dofs):

        cdef vector[int] global_dofs_vector
        cdef vector[vector[c_LocalDof]] local_dofs_vector
        cdef vector[vector[double]] weights_vector

        cdef int i
        cdef int j
        global_dofs_vector.resize(len(global_dofs))

        for i in range(len(global_dofs)):
            global_dofs_vector[i] = global_dofs[i]

        deref(self.impl_).global2localDofs(global_dofs_vector,
                local_dofs_vector, weights_vector)

        local_dofs = []
        weights = []

        for i in range(local_dofs_vector.size()):
            local_dofs.append([])
            for j in range(local_dofs_vector[i].size()):
                local_dofs[i].append(LocalDof(
                    local_dofs_vector[i][j].entity_index,
                    local_dofs_vector[i][j].dof_index))


        for i in range(weights_vector.size()):
            weights.append([])
            for j in range(weights_vector[i].size()):
                weights[i].append(weights_vector[i][j])

        return (local_dofs, weights)













    property global_dof_interpolation_points:
        """ (3xN) matrix of global interpolation points for the space, where each column is the
            coordinate of an interpolation points. """

        def __get__(self):
            cdef Matrix[double] data
            deref(self.impl_).getGlobalDofInterpolationPoints(data)
            return eigen_matrix_to_np_float64(data)

    property global_dof_normals:
        """ (3xN) matrix of normal directions associated with the interpolation points. """

        def __get__(self):
            cdef Matrix[double] data
            deref(self.impl_).getNormalsAtGlobalDofInterpolationPoints(data)
            return eigen_matrix_to_np_float64(data)

def function_space(Grid grid, kind, order, domains=None, cbool closed=True, cbool strictly_on_segment=False,
        cbool reference_point_on_segment=True, cbool element_on_segment=False,
        faces_to_include=None, nodes_to_include=None):
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
            "RT": Raviart-Thomas Vector spaces.
            "RWG": RWG Vector spaces.
            "NC": Nedelec Vector spaces.
            "SNC": Scaled Nedelec Vector spaces.

            "B-P": Polynomial spaces on barycentric grids.
            "B-DP": Polynomial discontinuous spaces on barycentric grids.
            "B-RT": Raviart-Thomas Vector spaces on barycentric grids.
            "B-RWG": RWG Vector spaces on barycentric grids.
            "B-NC": Nedelec Vector spaces on barycentric grids.
            "B-SNC": Scaled Nedelec Vector spaces on barycentric grids.

            "DUAL": Dual space on dual grid (only implemented for constants and linears).
            "BC": Buffa-Christian Vector space.
            "RBC": Rotated Buffa-Christiansen Vector space.

            "DP-CUSTOM" : Used internally

    order : int
        The order of the space, e.g. 0 for piecewise const, 1 for
        piecewise linear functions.

    domains : list
        List of integers specifying a list of physical entities
        of subdomains that should be included in the space.

    closed : bool
        Specifies whether the space is defined on a closed
        or open subspace.

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
    cdef Space s = Space()
    cdef Grid bary_grid
    cdef int dof_mode = 0

    if element_on_segment:
        dof_mode += 2
    if reference_point_on_segment:
        dof_mode += 4

    if dof_mode < 2:
        raise ValueError("At least one of 'reference_point_on_segment' or 'element_on_segment' must be true.")

    if kind=="P":
        if not (order>=1 and order <=10):
            raise ValueError("Order must be between 1 and 10")
        if (order==1):
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpace[double](grid.impl_))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
        else:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialContinuousScalarSpace[double](grid.impl_, order))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialContinuousScalarSpace[double](grid.impl_, order, domains, closed, strictly_on_segment))))
    elif kind=="DP":
        if not (order>=0 and order <=10):
            raise ValueError("Order must be between 0 and 10")
        if order==0:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpace[double](grid.impl_))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpace[double](grid.impl_, domains, closed))))
        else:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialDiscontinuousScalarSpace[double](grid.impl_, order))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialDiscontinuousScalarSpace[double](grid.impl_, order, domains, closed, dof_mode))))
    elif kind=="DP-CUSTOM":
        if order != 0 and order != 1:
            raise ValueError("Order must be 0 or 1")
        if order==0:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpaceOnCustomSegment[double](grid.impl_, faces_to_include, nodes_to_include))))
        if order == 1:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseLinearSpaceOnCustomSegment[double](grid.impl_, faces_to_include, nodes_to_include))))
    elif kind=="RT":
        if order!=0:
            raise ValueError("Only 0 order Raviart-Thomas spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="NC":
        if order!=0:
            raise ValueError("Only 0 order Nedelec spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveNedelec0VectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveNedelec0VectorSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="SNC":
        if order!=0:
            raise ValueError("Only 0 order Nedelec spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveScaledNedelec0VectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveScaledNedelec0VectorSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="RWG":
        if order!=0:
            raise ValueError("Only 0 order RWG spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRWGVectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRWGVectorSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="DUAL":
        if order == 0:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseConstantDualGridScalarSpace[double](grid.impl_))
                ))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseConstantDualGridScalarSpace[double](grid.impl_, domains, closed, strictly_on_segment))
                ))
        elif order == 1:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseLinearDualGridContinuousScalarSpace[double](grid.impl_))
                ))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseLinearDualGridContinuousScalarSpace[double](grid.impl_, domains, closed, strictly_on_segment))
                ))
        else:
            raise ValueError("Only order 0 and 1 dual grid spaces are implemented.")

    elif kind == "B-P":
        if order != 1:
            raise ValueError("Only linear spaces on barycentric grids are supported.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpaceBarycentric[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpaceBarycentric[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind == "B-DP":
        if order == 1:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseLinearDiscontinuousScalarSpaceBarycentric[double](grid.impl_))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseLinearDiscontinuousScalarSpaceBarycentric[double](grid.impl_, domains))))
        elif order == 0:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpaceBarycentric[double](grid.impl_))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptivePiecewiseConstantScalarSpaceBarycentric[double](grid.impl_, domains))))
        else:
            raise ValueError("Only constant and linear spaces on barycentric grids are supported.")

    elif kind=="B-RT":
        if order!=0:
            raise ValueError("Only 0 order Raviart-Thomas spaces on barycentric grids are supported.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpaceBarycentric[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpaceBarycentric[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="B-NC":
        if order!=0:
            raise ValueError("Only 0 order Nedelec spaces on barycentric grids are supported.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveNedelec0VectorSpaceBarycentric[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveNedelec0VectorSpaceBarycentric[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="B-SNC":
        if order!=0:
            raise ValueError("Only 0 order Nedelec spaces on barycentric grids are supported.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveScaledNedelec0VectorSpaceBarycentric[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveScaledNedelec0VectorSpaceBarycentric[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="B-RWG":
        if order!=0:
            raise ValueError("Only 0 order RWG spaces on barycentric grids are supported.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRWGVectorSpaceBarycentric[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRWGVectorSpaceBarycentric[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="BC":
        if order!=0:
            raise ValueError("Only order 0 Buffa-Christiansen spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveBuffaChristiansenVectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveBuffaChristiansenVectorSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
    elif kind=="RBC":
        if order!=0:
            raise ValueError("Only order 0 Buffa-Christiansen spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRotatedBuffaChristiansenVectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRotatedBuffaChristiansenVectorSpace[double](grid.impl_, domains, closed, strictly_on_segment))))
    else:
        raise ValueError("Unknown kind")

    return s

def evaluate_local_surface_gradient_ext(Space space, Entity0 element, object local_coordinates, object local_coefficients):
    """Evaluate the surface gradient on a given element."""

    import numpy as np

    if np.isreal(local_coefficients).all():
        coeffs_real = local_coefficients
        return eigen_matrix_to_np_float64(
                c_evaluateSurfaceGradients[double](deref(space.impl_), deref(element.impl_),
                                            np_to_eigen_matrix_float64(local_coordinates),
                                            np_to_eigen_vector_float64(local_coefficients)))
    else:
        coeffs_complex = local_coefficients
        return eigen_matrix_to_np_complex128(
                c_evaluateSurfaceGradients[complex_double](deref(space.impl_), deref(element.impl_),
                                                    np_to_eigen_matrix_float64(local_coordinates),
                                                    np_to_eigen_vector_complex128(local_coefficients)))
