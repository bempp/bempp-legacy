from unittest import TestCase


class TestBlockedBoundaryOperator(TestCase):

    def setUp(self):
        import bempp.api

        grid = bempp.api.shapes.regular_sphere(2)
        self._const_space = bempp.api.function_space(grid, "DP", 0)
        self._lin_space = bempp.api.function_space(grid, "P", 1)
        
        self._slp = bempp.api.operators.boundary.laplace.single_layer(self._const_space,
                self._lin_space, self._const_space)

        self._hyp = bempp.api.operators.boundary.laplace.hypersingular(self._lin_space,
                self._const_space, self._lin_space)

    def test_assembly_of_blocked_operator_fails_for_incompatible_dimensions(self):

        from bempp.api import BlockedOperator

        op = BlockedOperator(2,1)
        op[0, 0] = self._slp
        with self.assertRaises(ValueError):
            op[1, 0] = self._hyp

    def test_assembly_of_blocked_operator_succeeds_for_ompatible_dimensions(self):

        from bempp.api import BlockedOperator

        op = BlockedOperator(2,1)
        op[0, 0] = self._slp
        op[1, 0] = self._slp

        weak_op = op.weak_form()
        shape = weak_op.shape

        self.assertEqual(shape[0], 2 * self._slp.weak_form().shape[0])
        self.assertEqual(shape[1], self._slp.weak_form().shape[1])

    def test_matvec_with_pair_of_spaces(self):

        from bempp.api import BlockedOperator
        from bempp.api import GridFunction
        import bempp.api
        import numpy as np

        op = BlockedOperator(2, 2)
        for i in range(2):
            for j in range(2):
                op[i, j] = self._slp

        dof_count = self._const_space.global_dof_count
        grid_fun = bempp.api.GridFunction(self._const_space,
                coefficients=np.ones(dof_count))
        res = op * [grid_fun, grid_fun]
        proj = res[0].projections()

        single_res = self._slp * grid_fun

        np.testing.assert_allclose(res[0].projections(),
                2 * single_res.projections())










if __name__ == "__main__":
    from unittest import main

    main()
