"""Test the assembler routines."""

from unittest import TestCase

#pylint: disable=invalid-name

class TestAssembler(TestCase):
    """Test the assembler."""

    def setUp(self):
        """Setup test cases."""
        import bempp.api

        grid = bempp.api.shapes.regular_sphere(3)
        space = bempp.api.function_space(grid, "P", 1)
        pc_space = bempp.api.function_space(grid, "DP", 0)

        self._real_operator = \
            bempp.api.operators.boundary.laplace.single_layer(
                space, space, space, use_projection_spaces=False)
        self._real_operator_2 = \
            bempp.api.operators.boundary.laplace.single_layer(
                space, space, pc_space, use_projection_spaces=False)
        self._complex_operator = \
            bempp.api.operators.boundary.helmholtz.single_layer(
                space, space, space, 1, use_projection_spaces=False)

        self._rows = (5, 10)
        self._cols = (73, 100)

    def tearDown(self):
        """Finalise test cases."""
        import bempp.api
        bempp.api.global_parameters = bempp.api.common.global_parameters()

    def test_assemble_complete_dense_real_operator(self):
        """Assemble a dense real operator."""
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = \
            'dense'

        operator = self._real_operator

        actual = as_matrix(assemble_dense_block(
            self._real_operator, bempp.api.ALL, bempp.api.ALL))


        expected = as_matrix(operator.weak_form())

        self.assertEqual(258, actual.shape[0])
        self.assertEqual(258, actual.shape[1])
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_assemble_complete_dense_real_operator_non_square(self):
        "Assemble an opperator with a non square matrix"
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = \
            'dense'

        operator = self._real_operator_2

        actual = as_matrix(
            assemble_dense_block(operator, bempp.api.ALL, bempp.api.ALL))

        expected = as_matrix(operator.weak_form())

        self.assertEqual(512, actual.shape[0])
        self.assertEqual(258, actual.shape[1])
        np.testing.assert_allclose(actual, expected)

    def test_assemble_subblock_dense_real_operator(self):
        """Assemble a subblock of a dense real operator."""
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = \
            'dense'

        operator = self._real_operator

        actual = as_matrix(
            assemble_dense_block(self._real_operator, self._rows, self._cols))

        expected = as_matrix(operator.weak_form())[self._rows[
            0]:self._rows[1], self._cols[0]:self._cols[1]]

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_assemble_complete_dense_complex_operator(self):
        """Assemble a complex dense operator."""
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp.api

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = \
            'dense'

        operator = self._complex_operator
        actual = as_matrix(assemble_dense_block(
            self._complex_operator, bempp.api.ALL, bempp.api.ALL))

        expected = as_matrix(operator.weak_form())

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_assemble_subblock_dense_complex_operator(self):
        """Assembly a subblock of a complex dense operator."""
        from bempp.api import as_matrix, assemble_dense_block
        import numpy as np
        import bempp.api

        bempp.api.global_parameters.assembly.boundary_operator_assembly_type = \
            'dense'

        operator = self._complex_operator

        actual = as_matrix(
            assemble_dense_block(
                self._complex_operator, self._rows, self._cols))

        expected = as_matrix(operator.weak_form())[self._rows[
            0]:self._rows[1], self._cols[0]:self._cols[1]]

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    #pylint: disable=too-many-locals
    #pylint: disable=no-self-use
    def test_singular_assembly_of_laplace_single_layer_operator(self):
        """Assemble the singular part of the Laplace single-layer operator."""
        import bempp.api
        import numpy as np

        grid = bempp.api.shapes.cube(h=0.1)
        space = bempp.api.function_space(grid, "P", 1)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        op_sing = bempp.api.operators.boundary.laplace.single_layer(
            space, space, space,
            assemble_only_singular_part=True, use_projection_spaces=False)

        op = bempp.api.operators.boundary.laplace.single_layer(
            space, space, space,
            use_projection_spaces=False, parameters=parameters)

        mat_sing = op_sing.weak_form().sparse_operator

        mat = bempp.api.as_matrix(op.weak_form())

        indices = mat_sing.nonzero()
        n_elements = mat_sing.nnz

        actual = np.zeros(n_elements, dtype='float64')
        expected = np.zeros(n_elements, dtype='float64')

        for i in range(n_elements):

            m = indices[0][i]
            n = indices[1][i]

            actual[i] = mat_sing[m, n]
            expected[i] = mat[m, n]

        from numpy.testing import assert_allclose
        assert_allclose(actual, expected, rtol=1E-13)


    def test_singular_assembly_of_maxwell_efie_operator(self):
        """Assemble singular part of the Maxwell EFIE operator."""
        import bempp.api
        import numpy as np

        grid = bempp.api.shapes.cube(h=0.1)
        rt_space = bempp.api.function_space(grid, "RT", 0)
        nc_space = bempp.api.function_space(grid, "NC", 0)

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

        op_sing = bempp.api.operators.boundary.maxwell.electric_field(
            rt_space, rt_space, nc_space, 1,
            assemble_only_singular_part=True, use_projection_spaces=False)

        op = bempp.api.operators.boundary.maxwell.electric_field(
            rt_space, rt_space, nc_space, 1,
            use_projection_spaces=False, parameters=parameters)

        mat_sing = op_sing.weak_form().sparse_operator

        mat = bempp.api.as_matrix(op.weak_form())

        indices = mat_sing.nonzero()
        n_elements = mat_sing.nnz

        actual = np.zeros(n_elements, dtype='complex128')
        expected = np.zeros(n_elements, dtype='complex128')

        for i in range(n_elements):

            m = indices[0][i]
            n = indices[1][i]

            actual[i] = mat_sing[m, n]
            expected[i] = mat[m, n]

        self.assertAlmostEqual(np.linalg.norm(
            actual-expected, np.inf) / np.linalg.norm(expected, np.inf), 0)


if __name__ == "__main__":
    from unittest import main

    main()
