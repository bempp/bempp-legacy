from unittest import TestCase
import bempp

WAVE_NUMBER = -2j

class TestModifiedHelmholtz(TestCase):
    """Test cases for modified Helmholtz operators."""

    def setUp(self):
        grid = bempp.api.shapes.regular_sphere(3)
        self._lin_space = bempp.api.function_space(grid, "P", 1)
        self._const_space = bempp.api.function_space(grid, "DP", 0)

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


if __name__ == "__main__":
    from unittest import main

    main()
