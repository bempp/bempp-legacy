from unittest import TestCase
import bempp

WAVE_NUMBER = 1


class TestMaxwell(TestCase):
    """Test cases for Maxwell operators."""

    def setUp(self):
        grid = bempp.api.shapes.regular_sphere(3)
        self._space = bempp.api.function_space(grid, "RT", 0)

    def test_compound_electric_field_operator_agrees_with_standard_electric_field_operator(self):
        from bempp.api.operators.boundary.maxwell import electric_field
        import numpy as np

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        standard_efie = electric_field(self._space,
                                       WAVE_NUMBER,
                                       parameters=parameters).weak_form()
        compound_efie = electric_field(self._space,
                                       WAVE_NUMBER,
                                       parameters=parameters, use_slp=True).weak_form()

        mat1 = bempp.api.as_matrix(standard_efie)
        mat2 = bempp.api.as_matrix(compound_efie)

        self.assertAlmostEqual(np.linalg.norm(mat1 - mat2) / np.linalg.norm(mat1), 0)


if __name__ == "__main__":
    from unittest import main

    main()
