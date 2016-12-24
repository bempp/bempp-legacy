"""Test the modified Helmholtz operators."""

from unittest import TestCase
import bempp.api
import numpy as np

WAVE_NUMBER = -2j


class TestModifiedHelmholtz(TestCase):
    """Test cases for modified Helmholtz operators."""

    def setUp(self):
        """Setup unit tests."""
        self._grid = bempp.api.shapes.regular_sphere(2)
        self._lin_space = bempp.api.function_space(self._grid, "P", 1)
        self._const_space = bempp.api.function_space(self._grid, "DP", 0)

    def test_compound_hyp_agrees_with_standard_hyp_lin(self):
        """Compount HYP agrees with standard HYP for linear spaces."""
        from bempp.api.operators.boundary.modified_helmholtz import \
            hypersingular

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        standard_hypersingular = hypersingular(
            self._lin_space, self._lin_space, self._lin_space,
            WAVE_NUMBER,
            parameters=parameters, use_slp=False).weak_form()
        compound_hypersingular = hypersingular(
            self._lin_space, self._lin_space, self._lin_space,
            WAVE_NUMBER,
            parameters=parameters, use_slp=True).weak_form()

        mat1 = bempp.api.as_matrix(standard_hypersingular)
        mat2 = bempp.api.as_matrix(compound_hypersingular)

        self.assertAlmostEqual(np.linalg.norm(
            mat1 - mat2) / np.linalg.norm(mat1), 0)

    def test_multitrace_single_layer_agrees_with_standard_single_layer_dual_lin(
            self):
        """Multitrace SLP agrees with standard SLP."""
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        from bempp.api.operators.boundary.modified_helmholtz import \
            multitrace_operator, single_layer

        multitrace = multitrace_operator(
            self._grid, WAVE_NUMBER, parameters=parameters)

        actual = multitrace[0, 1]

        expected = single_layer(
            self._lin_space, self._lin_space, self._lin_space, WAVE_NUMBER,
            parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_multitrace_hypersingular_agrees_with_standard_hypersingular_lin(
            self):
        """Multitrace hyp agrees with standard hyp on linear spaces."""
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        from bempp.api.operators.boundary.modified_helmholtz import \
            multitrace_operator, hypersingular


        multitrace = multitrace_operator(
            self._grid, WAVE_NUMBER, parameters=parameters)

        actual = multitrace[1, 0]

        expected = hypersingular(
            self._lin_space, self._lin_space, self._lin_space, WAVE_NUMBER,
            parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 3)

    def test_multitrace_dlp_agrees_with_standard_dlp_lin(self):
        """Multitrace dlp agrees with standard dlp on linear spaces."""

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        from bempp.api.operators.boundary.modified_helmholtz import \
            multitrace_operator, double_layer


        multitrace = multitrace_operator(
            self._grid, WAVE_NUMBER, parameters=parameters)
        dlp = double_layer(
            self._lin_space, self._lin_space, self._lin_space, WAVE_NUMBER,
            parameters=parameters)

        actual = multitrace[0, 0]
        expected = -dlp

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_multitrace_adjoint_dlp_agrees_with_standard_adjoint_dlp_lin(self):
        """Multitrace ADLP agrees with standard ADLP on linear spaces."""

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        from bempp.api.operators.boundary.modified_helmholtz import \
            multitrace_operator, adjoint_double_layer


        multitrace = multitrace_operator(
            self._grid, WAVE_NUMBER, parameters=parameters)
        adlp = adjoint_double_layer(
            self._lin_space, self._lin_space, self._lin_space,
            WAVE_NUMBER, parameters=parameters)

        actual = multitrace[1, 1]
        expected = adlp

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 4)

    def test_compound_hypersingular_agrees_with_standard_hypersingular_operator(self):
        from bempp.api.operators.boundary.modified_helmholtz import hypersingular
        import numpy as np

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        standard_hypersingular = hypersingular(self._lin_space, self._const_space, self._lin_space,
                                               WAVE_NUMBER,
                                               parameters=parameters, use_slp=False).weak_form()
        compound_hypersingular = hypersingular(self._lin_space, self._const_space, self._lin_space,
                                               WAVE_NUMBER,
                                               parameters=parameters, use_slp=True).weak_form()

        mat1 = bempp.api.as_matrix(standard_hypersingular)
        mat2 = bempp.api.as_matrix(compound_hypersingular)

        self.assertAlmostEqual(np.linalg.norm(
            mat1 - mat2) / np.linalg.norm(mat1), 0)

    def test_multitrace_single_layer_agrees_with_standard_single_layer_dual(self):

        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        multitrace = bempp.api.operators.boundary.modified_helmholtz.multitrace_operator(self._grid, WAVE_NUMBER,
                                                                                         parameters=parameters,
                                                                                         spaces='dual')

        actual = multitrace[0, 1]

        expected = bempp.api.operators.boundary.modified_helmholtz.single_layer(dual_const_space, dual_const_space,
                                                                                dual_const_space, WAVE_NUMBER,
                                                                                parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_multitrace_hypersingular_agrees_with_standard_hypersingular_dual(self):

        lin_space = bempp.api.function_space(self._grid, "P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        multitrace = bempp.api.operators.boundary.modified_helmholtz.multitrace_operator(self._grid, WAVE_NUMBER,
                                                                                         parameters=parameters,
                                                                                         spaces='dual')

        actual = multitrace[1, 0]

        expected = bempp.api.operators.boundary.modified_helmholtz.hypersingular(lin_space, lin_space,
                                                                                 lin_space, WAVE_NUMBER,
                                                                                 parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 3)

    def test_multitrace_dlp_agrees_with_standard_dlp_dual(self):

        const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_bary = bempp.api.function_space(self._grid, "B-P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        multitrace = bempp.api.operators.boundary.modified_helmholtz.multitrace_operator(self._grid, WAVE_NUMBER,
                                                                                         parameters=parameters,
                                                                                         spaces='dual')
        dlp = bempp.api.operators.boundary.modified_helmholtz.double_layer(lin_space_bary, lin_space_bary, const_space,
                                                                           WAVE_NUMBER, parameters=parameters)
        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space_bary, lin_space_bary, const_space)

        actual = multitrace[0, 0]
        expected = -dlp

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_multitrace_adjoint_dlp_agrees_with_standard_adjoint_dlp_dual(self):

        const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_bary = bempp.api.function_space(self._grid, "B-P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        multitrace = bempp.api.operators.boundary.modified_helmholtz.multitrace_operator(self._grid, WAVE_NUMBER,
                                                                                         parameters=parameters,
                                                                                         spaces='dual')
        adlp = bempp.api.operators.boundary.modified_helmholtz.adjoint_double_layer(const_space, const_space, lin_space_bary,
                                                                                    WAVE_NUMBER, parameters=parameters)
        ident = bempp.api.operators.boundary.sparse.identity(
            const_space, const_space, lin_space_bary)

        actual = multitrace[1, 1]
        expected = adlp

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(
            mat_expected - mat_actual) / np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 4)

    def test_slp_hyp_pair_on_dual_grids(self):

        import bempp.api
        lin_space = self._lin_space
        const_space = self._const_space

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        expected_ops_dual = bempp.api.operators.boundary.modified_helmholtz.single_layer_and_hypersingular_pair(
            self._grid, WAVE_NUMBER, spaces='dual', parameters=parameters)

        dual_space = bempp.api.function_space(self._grid, "DUAL", 0)
        expected_slp = bempp.api.as_matrix(bempp.api.operators.boundary.modified_helmholtz.single_layer(
            dual_space, dual_space, dual_space, WAVE_NUMBER, parameters=parameters).weak_form())
        actual_slp = bempp.api.as_matrix(expected_ops_dual[0].weak_form())

        expected_hyp = bempp.api.as_matrix(bempp.api.operators.boundary.modified_helmholtz.hypersingular(
            lin_space, lin_space, lin_space, WAVE_NUMBER, parameters=parameters).weak_form())
        actual_hyp = bempp.api.as_matrix(expected_ops_dual[1].weak_form())

        diff_norm_slp = np.linalg.norm(
            expected_slp - actual_slp) / np.linalg.norm(actual_slp)
        diff_norm_hyp = np.linalg.norm(
            expected_hyp - actual_hyp) / np.linalg.norm(actual_hyp)

        self.assertAlmostEqual(diff_norm_slp, 0, 6)
        self.assertAlmostEqual(diff_norm_hyp, 0, 4)

    def test_slp_hyp_pair_linear_space(self):

        import bempp.api
        lin_space = self._lin_space
        const_space = self._const_space

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.quadrature.double_singular = 9
        parameters.quadrature.near.double_order = 9
        parameters.quadrature.medium.double_order = 9
        parameters.quadrature.far.double_order = 9

        expected_ops_dual = bempp.api.operators.boundary.modified_helmholtz.single_layer_and_hypersingular_pair(
            self._grid, WAVE_NUMBER, spaces='linear', parameters=parameters)

        expected_slp = bempp.api.as_matrix(bempp.api.operators.boundary.modified_helmholtz.single_layer(
            lin_space, lin_space, lin_space, WAVE_NUMBER, parameters=parameters).weak_form())
        actual_slp = bempp.api.as_matrix(expected_ops_dual[0].weak_form())

        expected_hyp = bempp.api.as_matrix(bempp.api.operators.boundary.modified_helmholtz.hypersingular(
            lin_space, lin_space, lin_space, WAVE_NUMBER, parameters=parameters).weak_form())
        actual_hyp = bempp.api.as_matrix(expected_ops_dual[1].weak_form())

        diff_norm_slp = np.linalg.norm(
            expected_slp - actual_slp) / np.linalg.norm(actual_slp)
        diff_norm_hyp = np.linalg.norm(
            expected_hyp - actual_hyp) / np.linalg.norm(actual_hyp)

        self.assertAlmostEqual(diff_norm_slp, 0, 6)
        self.assertAlmostEqual(diff_norm_hyp, 0, 6)


if __name__ == "__main__":
    from unittest import main

    main()
