"""Test P1 coupling with FEniCS."""

from unittest import TestCase
from unittest import skipIf
import bempp.api
import numpy as np

#pylint: disable=invalid-name

@skipIf(not bempp.api.HAVE_DOLFIN, "Test requires Dolfin")
class TestP1Coupling(TestCase):
    """Test coupling with FEniCS."""

    #pylint: disable=too-many-locals
    def test_dolfin_p1_identity_equals_bempp_p1_identity(self):
        """Dolfin P1 boundary identity is equal to BEM++ P1 identity."""

        #pylint: disable=import-error
        import dolfin
        from bempp.api.fenics_interface import fenics_to_bempp_trace_data
        from bempp.api.fenics_interface import FenicsOperator

        mesh = dolfin.UnitCubeMesh(5, 5, 5)
        V = dolfin.FunctionSpace(mesh, "CG", 1)

        space, trace_matrix = fenics_to_bempp_trace_data(V)

        u = dolfin.TestFunction(V)
        v = dolfin.TrialFunction(V)
        a = dolfin.inner(u, v) * dolfin.ds
        fenics_mass = FenicsOperator(a).weak_form().sparse_operator
        #pylint: disable=no-member
        actual = trace_matrix * fenics_mass * trace_matrix.transpose()

        from bempp.api.operators.boundary import sparse

        expected = sparse.identity(
            space, space, space).weak_form().sparse_operator
        diff = actual - expected
        self.assertAlmostEqual(np.max(np.abs(diff.data)), 0)


if __name__ == "__main__":
    from unittest import main

    main()
