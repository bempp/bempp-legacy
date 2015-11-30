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

def function_space(Grid grid, kind, order, domains=None, cbool closed=True):

    cdef Space s = Space()
    cdef Grid bary_grid
    if kind=="P":
        if not (order>=1 and order <=10):
            raise ValueError("Order must be between 1 and 10")
        if (order==1):
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpace[double](grid.impl_))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpace[double](grid.impl_, domains, closed))))
        else:
            if domains is None:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialContinuousScalarSpace[double](grid.impl_, order))))
            else:
                s.impl_.assign(reverse_const_pointer_cast(
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialContinuousScalarSpace[double](grid.impl_, order, domains, closed))))
    elif kind=="DP":
        if not (order>=0 and order <=10):
            raise ValueError("Order must be between 0 and 10")
        if (order==0):
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
                        shared_ptr[c_Space[double]](adaptivePiecewisePolynomialDiscontinuousScalarSpace[double](grid.impl_, order, domains, closed))))
    elif kind=="RT":
        if order!=0:
            raise ValueError("Only 0 order Raviart-Thomas spaces are implemented.")
        if domains is None:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpace[double](grid.impl_))))
        else:
            s.impl_.assign(reverse_const_pointer_cast(
                    shared_ptr[c_Space[double]](adaptiveRaviartThomas0VectorSpace[double](grid.impl_, domains, closed))))
    elif kind=="DUAL":
        if order != 0:
            raise ValueError("Only 0 order dual grid spaces are implemented.")
        if domains is not None:
            raise ValueError("Spaces on subdomains are not supported on dual grids.")
        s.impl_.assign(reverse_const_pointer_cast(
            shared_ptr[c_Space[double]](adaptivePiecewiseConstantDualGridScalarSpace[double](grid.impl_))
        ))
    elif kind == "B-P":
        if order != 1:
            raise ValueError("Only linear spaces on barycentric grids are supported.")
        if domains is not None:
            raise ValueError("Spaces on subdomains are not supported on barycentric grids.")
        s.impl_.assign(reverse_const_pointer_cast(
            shared_ptr[c_Space[double]](adaptivePiecewiseLinearContinuousScalarSpaceBarycentric[double](grid.impl_))))
    elif kind == "B-DP":
        if order != 1:
            raise ValueError("Only linear spaces on barycentric grids are supported.")
        if domains is not None:
            raise ValueError("Spaces on subdomains are not supported on barycentric grids.")
        s.impl_.assign(reverse_const_pointer_cast(
            shared_ptr[c_Space[double]](adaptivePiecewiseLinearDiscontinuousScalarSpaceBarycentric[double](grid.impl_))))
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


