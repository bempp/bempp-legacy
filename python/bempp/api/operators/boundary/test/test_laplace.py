from unittest import TestCase
from unittest import skip
import bempp.api
import numpy as np


class TestLaplace(TestCase):
    """Test cases for Laplace operators."""

    def setUp(self):
        grid = bempp.api.shapes.regular_sphere(2)
        self._lin_space = bempp.api.function_space(grid, "P", 1)
        self._const_space = bempp.api.function_space(grid, "DP", 0)
        self._grid = grid

    def test_compound_hypersingular_agrees_with_standard_hypersingular_operator(self):
        from bempp.api.operators.boundary.laplace import hypersingular
        import numpy as np

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        standard_hypersingular = hypersingular(self._lin_space, self._const_space, self._lin_space,
                                               parameters=parameters).weak_form()
        compound_hypersingular = hypersingular(self._lin_space, self._const_space, self._lin_space,
                                               parameters=parameters, use_slp=True).weak_form()

        mat1 = bempp.api.as_matrix(standard_hypersingular)
        mat2 = bempp.api.as_matrix(compound_hypersingular)

        self.assertAlmostEqual(np.linalg.norm(mat1 - mat2) / np.linalg.norm(mat1), 0)

    def test_dual_space_laplace_by_projection_is_correct(self):

        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_disc = bempp.api.function_space(self._grid.barycentric_grid(), "DP", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        slp = bempp.api.operators.boundary.laplace.single_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                                parameters=parameters)

        expected = bempp.api.operators.boundary.laplace.single_layer(dual_const_space, dual_const_space,
                                                                     dual_const_space,
                                                                     parameters=parameters)
        actual = bempp.api.project_operator(slp, domain=dual_const_space, range_=dual_const_space,
                                            dual_to_range=dual_const_space)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_linear_discontinuous_laplace_on_barycentric_by_projection_is_correct(self):

        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_disc = bempp.api.function_space(self._grid.barycentric_grid(), "DP", 1)
        lin_space_disc_bary = bempp.api.function_space(self._grid, "B-DP", 1)
        lin_space_disc_original = bempp.api.function_space(self._grid, "DP", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        slp = bempp.api.operators.boundary.laplace.single_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                                parameters=parameters)

        expected = bempp.api.operators.boundary.laplace.single_layer(lin_space_disc_original, lin_space_disc_original,
                                                                     lin_space_disc_original,
                                                                     parameters=parameters)
        actual = bempp.api.project_operator(slp, domain=lin_space_disc_bary, range_=lin_space_disc_bary,
                                            dual_to_range=lin_space_disc_bary)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 3)

    def test_hypersingular_on_barycentric_agrees_with_standard_hypersingular(self):

        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_disc = bempp.api.function_space(self._grid.barycentric_grid(), "DP", 1)
        lin_space_disc_bary = bempp.api.function_space(self._grid, "B-DP", 1)
        lin_space_bary = bempp.api.function_space(self._grid, "B-P", 1)
        lin_space = bempp.api.function_space(self._grid, "P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        slp = bempp.api.operators.boundary.laplace.single_layer(lin_space_disc, lin_space_disc, lin_space_disc,
                                                                parameters=parameters)

        slp_with_lin_disc_bary = bempp.api.project_operator(slp, domain=lin_space_disc_bary,
                                                            range_=lin_space_disc_bary,
                                                            dual_to_range=lin_space_disc_bary)

        expected = bempp.api.operators.boundary.laplace.hypersingular(lin_space, lin_space,
                                                                     lin_space,
                                                                     parameters=parameters)
        actual = bempp.api.operators.boundary.laplace.hypersingular(lin_space_bary, dual_const_space,
                                                                    lin_space_bary, use_slp=slp_with_lin_disc_bary,
                                                                    parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 4)


if __name__ == "__main__":
    from unittest import main

    main()
