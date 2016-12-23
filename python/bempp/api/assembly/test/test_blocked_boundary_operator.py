"""Test module for blocked boundary operator."""

from unittest import TestCase

#pylint: disable=invalid-name

class TestBlockedBoundaryOperator(TestCase):
    """Test blocked boundary operators."""

    def setUp(self):
        """Setup the test cases."""
        import bempp.api

        grid = bempp.api.shapes.regular_sphere(2)
        self._const_space = bempp.api.function_space(grid, "DP", 0)
        self._lin_space = bempp.api.function_space(grid, "P", 1)

        self._slp = bempp.api.operators.boundary.laplace.single_layer(
            self._const_space, self._lin_space, self._const_space)

        self._hyp = bempp.api.operators.boundary.laplace.hypersingular(
            self._lin_space, self._const_space, self._lin_space)

    def test_assembly_fails_for_incompatible_dimensions(self):
        """Assembly fails for incompatible dimensions."""

        from bempp.api import BlockedOperator

        op = BlockedOperator(2, 1)
        op[0, 0] = self._slp
        with self.assertRaises(ValueError):
            op[1, 0] = self._hyp

    def test_assembly_succeeds_for_compatible_dimensions(self):
        """Assembly succeeds for compatible dimensions."""

        from bempp.api import BlockedOperator

        op = BlockedOperator(2, 1)
        op[0, 0] = self._slp
        op[1, 0] = self._slp

        weak_op = op.weak_form()
        shape = weak_op.shape

        self.assertEqual(shape[0], 2 * self._slp.weak_form().shape[0])
        self.assertEqual(shape[1], self._slp.weak_form().shape[1])

    def test_matvec_with_pair_of_spaces(self):
        """Matvec gives the correct result."""

        from bempp.api import BlockedOperator
        import bempp.api
        import numpy as np

        op = BlockedOperator(2, 2)
        for i in range(2):
            for j in range(2):
                op[i, j] = self._slp

        dof_count = self._const_space.global_dof_count
        grid_fun = bempp.api.GridFunction(
            self._const_space, coefficients=np.ones(dof_count))
        res = op * [grid_fun, grid_fun]

        single_res = self._slp * grid_fun

        #pylint: disable=no-member
        np.testing.assert_allclose(
            res[0].projections(), 2 * single_res.projections())


if __name__ == "__main__":
    from unittest import main

    main()
