"""Wrapper for a FEniCS Operator."""

#pylint: disable=import-error
from dolfin import as_backend_type, assemble, parameters
from scipy.sparse import csr_matrix

#pylint: disable=too-few-public-methods
class FenicsOperator(object):
    """Wraps a FEniCS Operator into a BEM++ operator."""

    def __init__(self, fenics_weak_form):
        """Construct an operator from a weak form in FEniCS."""
        self._fenics_weak_form = fenics_weak_form
        self._sparse_mat = None

    def weak_form(self):
        """Return the weak form."""
        from bempp.api.assembly.discrete_boundary_operator import \
            SparseDiscreteBoundaryOperator

        if self._sparse_mat is not None:
            return SparseDiscreteBoundaryOperator(self._sparse_mat)

        backend = parameters['linear_algebra_backend']
        if backend not in ['PETSc', 'uBLAS']:
            parameters['linear_algebra_backend'] = 'PETSc'

        if parameters['linear_algebra_backend'] == 'PETSc':
            mat = as_backend_type(assemble(self._fenics_weak_form)).mat()
            (indptr, indices, data) = mat.getValuesCSR()
            self._sparse_mat = csr_matrix(
                (data, indices, indptr), shape=mat.size)
        elif parameters['linear_algebra_backend'] == 'uBLAS':
            self._sparse_mat = assemble(self._fenics_weak_form).sparray()
        else:
            raise ValueError(
                "This should not happen! Backend type is %s. " +
                "Default backend should have been set to 'PETSc'.",
                str(parameters['linear_algebra_backend']))

        parameters['linear_algebra_backend'] = backend
        return SparseDiscreteBoundaryOperator(self._sparse_mat)
