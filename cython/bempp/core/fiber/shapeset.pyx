from bempp.core.utils cimport np_to_eigen_matrix_float64
from bempp.core.fiber cimport _3d_array_to_numpy
from bempp.core.fiber cimport _4d_array_to_numpy


cdef class Shapeset:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __deallocate__(self):
        pass

    def evaluate(self, points, int dof, values=True, derivatives=False):
        """Evaluate the shapeset."""

        cdef Matrix[double] points_matrix = np_to_eigen_matrix_float64(
                points)

        cdef int what = 0

        cdef c_BasisData basis_data

        if values:
            what += 1
        if derivatives:
            what += 2

        if what == 0:
            return None

        self.impl_.evaluate(what, points_matrix, dof, basis_data)

        if what == 1:
            return _3d_array_to_numpy(basis_data.values)
        if what == 3:
            return (_3d_array_to_numpy(basis_data.values),
                    _4d_array_to_numpy(basis_data.derivatives))

    property order:

        def __get__(self):
            return self.impl_.order()

    property size:

        def __get__(self):
            return self.impl_.size()
