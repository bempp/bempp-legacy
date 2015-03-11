import pytest
import numpy as np

class TestP1Coupling(object):

    def test_dolfin_p1_identity_equals_bempp_p1_identity(self):

        import dolfin
        from bempp.fenics_interface import coupling

        mesh = dolfin.UnitCubeMesh(5,5,5)
        V = dolfin.FunctionSpace(mesh,"CG",1)

        space, trace_matrix = coupling.fenics_to_bempp_trace_data(V)

        u = dolfin.TestFunction(V)
        v = dolfin.TrialFunction(V)
        a = dolfin.inner(u,v)*dolfin.ds
        fenics_mass = dolfin.assemble(a).sparray()
        actual = trace_matrix*fenics_mass*trace_matrix.transpose()

        from bempp.operators.boundary import sparse

        expected = sparse.identity(space,space,space).weak_form().sparse_operator
        diff = actual-expected
        assert np.max(np.abs(diff.data))<1E-14
        


	
