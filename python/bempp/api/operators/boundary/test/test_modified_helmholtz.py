from unittest import TestCase
import bempp.api
import numpy as np

WAVE_NUMBER = -2j

class TestModifiedHelmholtz(TestCase):
    """Test cases for modified Helmholtz operators."""

    def setUp(self):
        self._grid = bempp.api.shapes.regular_sphere(2)
        self._lin_space = bempp.api.function_space(self._grid, "P", 1)
        self._const_space = bempp.api.function_space(self._grid, "DP", 0)

    def test_compound_hypersingular_agrees_with_standard_hypersingular_operator(self):
        from bempp.api.operators.boundary.modified_helmholtz import hypersingular
        import numpy as np

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        standard_hypersingular = hypersingular(self._lin_space, self._const_space, self._lin_space,
                                               WAVE_NUMBER,
                                               parameters=parameters).weak_form()
        compound_hypersingular = hypersingular(self._lin_space, self._const_space, self._lin_space,
                                               WAVE_NUMBER,
                                               parameters=parameters, use_slp=True).weak_form()

        mat1 = bempp.api.as_matrix(standard_hypersingular)
        mat2 = bempp.api.as_matrix(compound_hypersingular)

        self.assertAlmostEqual(np.linalg.norm(mat1 - mat2) / np.linalg.norm(mat1), 0)

    def test_calderon_single_layer_agrees_with_standard_single_layer(self):

        dual_const_space = bempp.api.function_space(self._grid, "DUAL", 0)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        calderon = bempp.api.operators.boundary.modified_helmholtz.interior_calderon_projector(self._grid, WAVE_NUMBER,
                                                                                               parameters=parameters)

        actual = calderon[0, 1]

        expected = bempp.api.operators.boundary.modified_helmholtz.single_layer(dual_const_space, dual_const_space,
                                                                     dual_const_space, WAVE_NUMBER,
                                                                     parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_calderon_hypersingular_agrees_with_standard_hypersingular(self):

        lin_space = bempp.api.function_space(self._grid, "P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        calderon = bempp.api.operators.boundary.modified_helmholtz.interior_calderon_projector(self._grid, WAVE_NUMBER,
                                                                                               parameters=parameters)

        actual = calderon[1, 0]

        expected = bempp.api.operators.boundary.modified_helmholtz.hypersingular(lin_space, lin_space,
                                                                      lin_space, WAVE_NUMBER,
                                                                      parameters=parameters)

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 3)

    def test_calderon_dlp_agrees_with_standard_dlp(self):

        const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_bary = bempp.api.function_space(self._grid, "B-P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        calderon = bempp.api.operators.boundary.modified_helmholtz.interior_calderon_projector(self._grid, WAVE_NUMBER,
                                                                                               parameters=parameters)
        dlp = bempp.api.operators.boundary.modified_helmholtz.double_layer(lin_space_bary, lin_space_bary, const_space,
                                                                           WAVE_NUMBER, parameters=parameters)
        ident = bempp.api.operators.boundary.sparse.identity(lin_space_bary, lin_space_bary, const_space)

        actual = calderon[0, 0]
        expected = .5 * ident - dlp

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0)

    def test_calderon_adjoint_dlp_agrees_with_standard_adjoint_dlp(self):

        const_space = bempp.api.function_space(self._grid, "DUAL", 0)
        lin_space_bary = bempp.api.function_space(self._grid, "B-P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        calderon = bempp.api.operators.boundary.modified_helmholtz.interior_calderon_projector(self._grid, WAVE_NUMBER,
                                                                                               parameters=parameters)
        adlp = bempp.api.operators.boundary.modified_helmholtz.adjoint_double_layer(const_space, const_space, lin_space_bary,
                                                                         WAVE_NUMBER, parameters=parameters)
        ident = bempp.api.operators.boundary.sparse.identity(const_space, const_space, lin_space_bary)

        actual = calderon[1, 1]
        expected = .5 * ident + adlp

        mat_expected = bempp.api.as_matrix(expected.weak_form())
        mat_actual = bempp.api.as_matrix(actual.weak_form())

        diff_norm = np.linalg.norm(mat_expected-mat_actual)/np.linalg.norm(mat_actual)

        self.assertAlmostEqual(diff_norm, 0, 4)


if __name__ == "__main__":
    from unittest import main

    main()
