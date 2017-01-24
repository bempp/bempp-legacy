"""Laplace operators unit tests."""

from unittest import TestCase
import bempp.api
import numpy as np


class TestLaplace(TestCase):
    """Test cases for Laplace operators."""

    def setUp(self):
        """Setup the unit tests."""
        grid = bempp.api.shapes.regular_sphere(2)
        self._lin_space = bempp.api.function_space(grid, "P", 1)
        self._const_space = bempp.api.function_space(grid, "DP", 0)
        self._grid = grid

    #pylint: disable=invalid-name
    def test_compound_hypersingular_agrees_with_standard_hypersingular_operator(
            self):
        """Compound and standard hyp agree."""
        from bempp.api.operators.boundary.laplace import hypersingular

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        standard_hypersingular = hypersingular(
            self._lin_space, self._const_space, self._lin_space,
            parameters=parameters, use_slp=False).weak_form()
        compound_hypersingular = hypersingular(
            self._lin_space, self._const_space, self._lin_space,
            parameters=parameters, use_slp=True).weak_form()

        mat1 = bempp.api.as_matrix(standard_hypersingular)
        mat2 = bempp.api.as_matrix(compound_hypersingular)

        self.assertAlmostEqual(np.linalg.norm(
            mat1 - mat2) / np.linalg.norm(mat1), 0)

    def test_dual_space_laplace_by_projection_is_correct(self):
        """Dual Space Laplace by projection is correct."""

        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_disc = bempp.api.function_space(
            self._grid.barycentric_grid(), "DP", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        slp = bempp.api.operators.boundary.laplace.single_layer(
            lin_space_disc, lin_space_disc, lin_space_disc,
            parameters=parameters)

        expected = bempp.api.operators.boundary.laplace.single_layer(
            dual_const_space, dual_const_space,
            dual_const_space,
            parameters=parameters)
        actual = bempp.api.project_operator(
            slp, domain=dual_const_space, range_=dual_const_space,
            dual_to_range=dual_const_space)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_linear_discontinuous_laplace_on_barycentric_by_proj_is_correct(
            self):
        """Linear disc. Laplace on barycentrics by projection is correct."""
        lin_space_disc = bempp.api.function_space(
            self._grid.barycentric_grid(), "DP", 1)
        lin_space_disc_bary = bempp.api.function_space(self._grid, "B-DP", 1)
        lin_space_disc_original = bempp.api.function_space(self._grid, "DP", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        slp = bempp.api.operators.boundary.laplace.single_layer(
            lin_space_disc, lin_space_disc, lin_space_disc,
            parameters=parameters)

        expected = bempp.api.operators.boundary.laplace.single_layer(
            lin_space_disc_original, lin_space_disc_original,
            lin_space_disc_original,
            parameters=parameters)
        actual = bempp.api.project_operator(
            slp, domain=lin_space_disc_bary, range_=lin_space_disc_bary,
            dual_to_range=lin_space_disc_bary)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 3)

    def test_hypersingular_on_barycentric_agrees_with_standard_hypersingular(
            self):
        """Hypersingular on barycentric agrees with standard hypersingular."""
        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_disc = bempp.api.function_space(
            self._grid.barycentric_grid(), "DP", 1)
        lin_space_disc_bary = bempp.api.function_space(self._grid, "B-DP", 1)
        lin_space_bary = bempp.api.function_space(self._grid, "B-P", 1)
        lin_space = bempp.api.function_space(self._grid, "P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        slp = bempp.api.operators.boundary.laplace.single_layer(
            lin_space_disc, lin_space_disc, lin_space_disc,
            parameters=parameters)

        slp_with_lin_disc_bary = bempp.api.project_operator(
            slp, domain=lin_space_disc_bary,
            range_=lin_space_disc_bary,
            dual_to_range=lin_space_disc_bary)

        expected = bempp.api.operators.boundary.laplace.hypersingular(
            lin_space, lin_space, lin_space,
            parameters=parameters)
        actual = bempp.api.operators.boundary.laplace.hypersingular(
            lin_space_bary, dual_const_space,
            lin_space_bary, use_slp=slp_with_lin_disc_bary,
            parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 4)

    def test_slp_hyp_pair_on_dual_grids(self):
        """SLP and HYP Pair on dual grids."""

        lin_space = self._lin_space

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        from bempp.api.operators.boundary.laplace import \
            single_layer_and_hypersingular_pair, single_layer, hypersingular

        actual_ops_dual = single_layer_and_hypersingular_pair(
            self._grid, spaces='dual', stabilization_factor=0,
            parameters=parameters)

        dual_space = bempp.api.function_space(self._grid, "DUAL", 0)
        expected_slp = bempp.api.as_matrix(
            single_layer(
                dual_space, dual_space, dual_space,
                parameters=parameters).weak_form())
        actual_slp = bempp.api.as_matrix(
            actual_ops_dual[0].weak_form())

        expected_hyp = bempp.api.as_matrix(
            hypersingular(
                lin_space, lin_space, lin_space,
                parameters=parameters).weak_form())
        actual_hyp = bempp.api.as_matrix(actual_ops_dual[1].weak_form())

        diff_norm_slp = np.linalg.norm(
            expected_slp - actual_slp) / np.linalg.norm(actual_slp)
        diff_norm_hyp = np.linalg.norm(
            expected_hyp - actual_hyp) / np.linalg.norm(actual_hyp)

        self.assertAlmostEqual(diff_norm_slp, 0, 6)
        self.assertAlmostEqual(diff_norm_hyp, 0, 4)

    def test_slp_hyp_pair_linear_space(self):
        """SLP and HYP Pair with linear spaces."""
        lin_space = self._lin_space

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        from bempp.api.operators.boundary.laplace import \
            single_layer_and_hypersingular_pair, single_layer, hypersingular

        actual_ops_dual = single_layer_and_hypersingular_pair(
            self._grid, spaces='linear', stabilization_factor=0,
            parameters=parameters)

        expected_slp = bempp.api.as_matrix(
            single_layer(
                lin_space, lin_space, lin_space,
                parameters=parameters).weak_form())
        actual_slp = bempp.api.as_matrix(actual_ops_dual[0].weak_form())

        expected_hyp = bempp.api.as_matrix(hypersingular(
            lin_space, lin_space, lin_space, parameters=parameters).weak_form())
        actual_hyp = bempp.api.as_matrix(actual_ops_dual[1].weak_form())

        diff_norm_slp = np.linalg.norm(
            expected_slp - actual_slp) / np.linalg.norm(actual_slp)
        diff_norm_hyp = np.linalg.norm(
            expected_hyp - actual_hyp) / np.linalg.norm(actual_hyp)

        self.assertAlmostEqual(diff_norm_slp, 0, 6)
        self.assertAlmostEqual(diff_norm_hyp, 0, 6)

    def test_regularized_hypersingular(self):
        """Regularized hypersingular operator."""
        lin_space = self._lin_space

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        from bempp.api.operators.boundary.laplace import \
            single_layer_and_hypersingular_pair

        actual_ops_dual = single_layer_and_hypersingular_pair(
            self._grid, spaces='linear', stabilization_factor=.5,
            parameters=parameters)

        expected_hyp = bempp.api.operators.boundary.laplace.hypersingular(
            lin_space, lin_space, lin_space, parameters=parameters)
        expected_hyp += 0.5 * \
            bempp.api.assembly.RankOneBoundaryOperator(
                lin_space, lin_space, lin_space)
        expected_hyp = bempp.api.as_matrix(expected_hyp.weak_form())
        actual_hyp = bempp.api.as_matrix(actual_ops_dual[1].weak_form())

        diff_norm_hyp = np.linalg.norm(
            expected_hyp - actual_hyp) / np.linalg.norm(actual_hyp)

        self.assertAlmostEqual(diff_norm_hyp, 0, 6)

if __name__ == "__main__":
    from unittest import main
    main()
