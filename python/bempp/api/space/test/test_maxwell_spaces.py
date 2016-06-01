from unittest import TestCase
import unittest
import bempp.api
import numpy as np


class TestMaxwellSpaces(TestCase):
    """Test Maxwell spaces."""

    requiresgmsh = unittest.skipIf(
        bempp.api.GMSH_PATH is None, reason="Needs GMSH")

    def setUp(self):

        self._grid = bempp.api.shapes.sphere(h=0.2)

    @requiresgmsh
    def test_rwg_space_is_correctly_scaled(self):

        space_rt = bempp.api.function_space(self._grid, "RT", 0)
        space_rwg = bempp.api.function_space(self._grid, "RWG", 0)

        for element in self._grid.leaf_view.entity_iterator(0):
            dofs_rt, weights_rt = space_rt.get_global_dofs(
                element, dof_weights=True)
            dofs_rwg, weights_rwg = space_rwg.get_global_dofs(
                element, dof_weights=True)
            self.assertEqual(len(weights_rt), len(weights_rwg))
            edge_volumes = np.array([edge.geometry.volume
                                     for edge in element.sub_entity_iterator(1)])
            self.assertEqual(dofs_rt, dofs_rwg)
            for ind in range(3):
                self.assertAlmostEqual(
                    weights_rt[ind] * edge_volumes[ind], weights_rwg[ind])

    @requiresgmsh
    def test_maxwell_barycentric_rwg_mass_agrees_with_rwg_mass(self):

        import numpy as np

        space_rwg = bempp.api.function_space(self._grid, "RWG", 0)
        space_rwg_b = bempp.api.function_space(self._grid, "B-RWG", 0)

        ident_rwg = bempp.api.operators.boundary.sparse.identity(
            space_rwg, space_rwg, space_rwg)
        ident_rwg_b = bempp.api.operators.boundary.sparse.identity(
            space_rwg_b, space_rwg_b, space_rwg_b)

        mat_ident_rwg = ident_rwg.weak_form().sparse_operator
        mat_ident_rwg_b = ident_rwg_b.weak_form().sparse_operator
        max_diff = np.max(np.abs((mat_ident_rwg - mat_ident_rwg_b).data))
        self.assertAlmostEqual(max_diff, 0)

if __name__ == "__main__":
    from unittest import main

    main()
