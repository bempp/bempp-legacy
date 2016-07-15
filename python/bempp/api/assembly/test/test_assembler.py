from unittest import TestCase


class TestAssembler(TestCase):

    def setUp(self):
        import bempp.api

        grid = bempp.api.shapes.regular_sphere(3)
        space = bempp.api.function_space(grid, "P", 1)
        pc_space = bempp.api.function_space(grid, "DP", 0)
        rt_space = bempp.api.function_space(grid, "RT", 0)
        nc_space = bempp.api.function_space(grid, "NC", 0)

        bempp.api.global_parameters.assembly.use_super_spaces = False

        self._real_operator = bempp.api.operators.boundary.laplace.single_layer(
            space, space, space)
        self._real_operator_2 = bempp.api.operators.boundary.laplace.single_layer(
            space, space, pc_space)
        self._complex_operator = bempp.api.operators.boundary.maxwell.electric_field(
            rt_space, rt_space, nc_space, 1)

        bempp.api.global_parameters.assembly.use_super_spaces = True

        self._rows = (5, 10)
        self._cols = (73, 100)

    def tearDown(self):
        import bempp

        bempp.api.global_parameters = bempp.api.common.global_parameters()

    def test_assemble_complete_dense_real_operator(self):
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'
        bempp.api.global_parameters.assembly.use_super_spaces = False

        operator = self._real_operator

        actual = as_matrix(assemble_dense_block(self._real_operator, bempp.api.ALL, bempp.api.ALL,
                                                operator.domain, operator.dual_to_range))

        bempp.api.global_parameters.assembly.use_super_spaces = True

        expected = as_matrix(operator.weak_form())

        self.assertEqual(258, actual.shape[0])
        self.assertEqual(258, actual.shape[1])
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_assemble_complete_dense_real_operator_non_square(self):
        "assemble an opperator with a non square matrix"
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'
        bempp.api.global_parameters.assembly.use_super_spaces = False

        operator = self._real_operator_2

        actual = as_matrix(assemble_dense_block(operator, bempp.api.ALL, bempp.api.ALL,
                                                operator.domain, operator.dual_to_range))
        bempp.api.global_parameters.assembly.use_super_spaces = True

        expected = as_matrix(operator.weak_form())

        self.assertEqual(512, actual.shape[0])
        self.assertEqual(258, actual.shape[1])
        np.testing.assert_allclose(actual, expected)

    def test_assemble_subblock_dense_real_operator(self):
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'

        operator = self._real_operator
        bempp.api.global_parameters.assembly.use_super_spaces = False

        actual = as_matrix(assemble_dense_block(self._real_operator, self._rows, self._cols,
                                                operator.domain, operator.dual_to_range))
        bempp.api.global_parameters.assembly.use_super_spaces = True

        expected = as_matrix(operator.weak_form())[self._rows[
            0]:self._rows[1], self._cols[0]:self._cols[1]]

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_assemble_complete_dense_complex_operator(self):
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'

        operator = self._complex_operator

        actual = as_matrix(assemble_dense_block(self._complex_operator, bempp.api.ALL, bempp.api.ALL,
                                                operator.domain, operator.dual_to_range))

        expected = as_matrix(operator.weak_form())

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_assemble_subblock_dense_complex_operator(self):
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'

        operator = self._complex_operator

        actual = as_matrix(assemble_dense_block(self._complex_operator, self._rows, self._cols,
                                                operator.domain, operator.dual_to_range))

        expected = as_matrix(operator.weak_form())[self._rows[
            0]:self._rows[1], self._cols[0]:self._cols[1]]

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

if __name__ == "__main__":
    from unittest import main

    main()
