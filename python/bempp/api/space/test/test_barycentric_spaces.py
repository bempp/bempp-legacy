from unittest import TestCase
import unittest
import bempp.api


class TestBarycentricLinearSpace(TestCase):
    """Test linear spaces on barycentric grids."""

    requiresgmsh = unittest.skipIf(bempp.api.GMSH_PATH is None, reason="Needs GMSH")

    def setUp(self):

        self._grid = bempp.api.shapes.sphere(h=0.2)
        self._bary_space = bempp.api.function_space(self._grid, "B-P", 1)
        self._space = bempp.api.function_space(self._grid, "P", 1)

    @requiresgmsh
    def test_global_dof_count_agrees_with_non_barycentric_space(self):

        self.assertEqual(self._bary_space.global_dof_count, self._space.global_dof_count)

    @requiresgmsh
    def test_mass_matrix_of_barycentric_space_agrees_with_non_barycentric_space(self):

        from bempp.api.operators.boundary.sparse import identity

        barycentric_ident = identity(self._bary_space, self._bary_space, self._bary_space).weak_form().sparse_operator
        ident = identity(self._space, self._space, self._space).weak_form().sparse_operator

        diff = barycentric_ident - ident

        self.assertAlmostEqual(diff.max(), 0, 15)
        self.assertAlmostEqual(diff.min(), 0, 15)

    @requiresgmsh
    def test_slp_operator_on_barycentric_space_agrees_with_slp_on_non_barycentric_space(self):

        import numpy as np
        slp = bempp.api.operators.boundary.laplace.single_layer
        space = self._space
        bary_space = self._bary_space

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        discrete_barycentric_slp = bempp.api.as_matrix(slp(space, space, space, parameters=parameters).weak_form())
        discrete_slp = bempp.api.as_matrix(slp(bary_space, bary_space, bary_space, parameters=parameters).weak_form())

        diff = np.linalg.norm(discrete_barycentric_slp - discrete_slp, np.inf) / np.linalg.norm(discrete_slp, np.inf)

        self.assertAlmostEqual(diff, 0, 4)

    @requiresgmsh
    def test_mass_matrix_of_barycentric_discontinuous_space_agrees_with_non_barycentric_discontinuous_space(self):

        from bempp.api.operators.boundary.sparse import identity

        disc_space_bary = bempp.api.function_space(self._grid, "B-DP", 1)
        disc_space = bempp.api.function_space(self._grid, "DP", 1)


        barycentric_ident = identity(disc_space_bary, disc_space_bary,
                                     disc_space_bary).weak_form().sparse_operator
        ident = identity(disc_space, disc_space,
                         disc_space).weak_form().sparse_operator

        diff = barycentric_ident - ident

        self.assertAlmostEqual(diff.max(), 0, 15)
        self.assertAlmostEqual(diff.min(), 0, 15)

    @requiresgmsh
    def test_slp_operator_on_barycentric_discontinuous_space_agrees_with_slp_on_non_barycentric_disc_space(self):

        import numpy as np
        slp = bempp.api.operators.boundary.laplace.single_layer
        bary_space = bempp.api.function_space(self._grid, "B-DP", 1)
        space = bempp.api.function_space(self._grid, "DP", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        discrete_barycentric_slp = bempp.api.as_matrix(slp(space, space, space, parameters=parameters).weak_form())
        discrete_slp = bempp.api.as_matrix(slp(bary_space, bary_space, bary_space, parameters=parameters).weak_form())

        diff = np.linalg.norm(discrete_barycentric_slp - discrete_slp, np.inf) / np.linalg.norm(discrete_slp, np.inf)

        self.assertAlmostEqual(diff, 0, 3)





if __name__ == "__main__":
    from unittest import main

    main()
