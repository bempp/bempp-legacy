from dolfin import as_backend_type, assemble, parameters
from scipy.sparse import csr_matrix


class FenicsOperator(object):

    def __init__(self, fenics_weak_form):
        self._fenics_weak_form = fenics_weak_form
        self._sparse_mat = None

    def weak_form(self):

        from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        if self._sparse_mat is not None:
            return SparseDiscreteBoundaryOperator(self._sparse_mat)

        backend = parameters['linear_algebra_backend']
        if backend not in ['PETSc', 'uBLAS']:
            parameters['linear_algebra_backend'] = 'PETSc'

        if parameters['linear_algebra_backend'] == 'PETSc':
            A = as_backend_type(assemble(self._fenics_weak_form)).mat()
            (indptr, indices, data) = A.getValuesCSR()
            self._sparse_mat = csr_matrix(
                (data, indices, indptr), shape=A.size)
        elif parameters['linear_algebra_backend'] == 'uBLAS':
            self._sparse_mat = assemble(self._fenics_weak_form).sparray()
        else:
            raise ValueError(
                "This should not happen! Backend type is '{0}'. " +
                "Default backend should have been set to 'PETSc'.".format(
                    parameters['linear_algebra_backend']))

        parameters['linear_algebra_backend'] = backend
        return SparseDiscreteBoundaryOperator(self._sparse_mat)
