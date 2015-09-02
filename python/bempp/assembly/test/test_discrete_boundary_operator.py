from unittest import TestCase


class TestGeneralNonlocalDiscreteBoundaryOperator(TestCase):
    """Test cases for general discrete operators."""

    def setUp(self):
        import bempp

        parameters = bempp.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'hmat'
        self._parameters = parameters
        grid = bempp.grid_from_sphere(3)
        space = bempp.function_space(grid, "DP", 0)
        self._operator_real = bempp.operators.boundary.laplace.single_layer(space, space, space,
                                                                            parameters=parameters).weak_form()
        self._operator_complex = bempp.operators.boundary.helmholtz.single_layer(space, space, space, 1,
                                                                                 parameters=parameters).weak_form()
        self._space = space

    def test_type(self):
        from bempp.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator

        self.assertIsInstance(self._operator_real, GeneralNonlocalDiscreteBoundaryOperator)

    def test_shape(self):
        expected = (self._space.global_dof_count, self._space.global_dof_count)
        actual = self._operator_real.shape
        self.assertEqual(expected, actual)

    def test_sum_of_real_and_complex_has_complex_type(self):
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator_real + self._operator_complex
        self.assertEqual(sum_operator.dtype, 'complex128')

    def matrix_vector_product_of_operator_sum_is_correct(self):
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator_real + self._operator_complex
        actual = sum_operator * vec
        expected = self._operator_real * vec + self._operator_complex * vec

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)


class TestDenseDiscreteBoundaryOperator(TestCase):
    """Test cases for dense discrete operators."""

    def setUp(self):
        import bempp

        parameters = bempp.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        self._parameters = parameters
        grid = bempp.grid_from_sphere(3)
        space = bempp.function_space(grid, "DP", 0)
        self._operator_real = bempp.operators.boundary.laplace.single_layer(space, space, space,
                                                                            parameters=parameters).weak_form()
        self._operator_complex = bempp.operators.boundary.helmholtz.single_layer(space, space, space, 1,
                                                                                 parameters=parameters).weak_form()
        self._space = space

    def test_type(self):
        from bempp.assembly.discrete_boundary_operator import DenseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator_real, DenseDiscreteBoundaryOperator)

    def test_shape(self):
        expected = (self._space.global_dof_count, self._space.global_dof_count)
        actual = self._operator_real.shape
        self.assertEqual(expected, actual)

    def test_sum_of_real_and_complex_has_complex_type(self):
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator_real + self._operator_complex
        self.assertEqual(sum_operator.dtype, 'complex128')

    def test_sum_of_dense_operators_is_dense(self):
        from bempp.assembly.discrete_boundary_operator import DenseDiscreteBoundaryOperator

        sum_operator = self._operator_real + self._operator_complex
        self.assertIsInstance(sum_operator, DenseDiscreteBoundaryOperator)

    def matrix_vector_product_of_operator_sum_is_correct(self):
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator_real + self._operator_complex
        actual = sum_operator * vec
        expected = self._operator_real * vec + self._operator_complex * vec

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_scalar_multiply_of_dense_operator_is_dense(self):
        from bempp.assembly.discrete_boundary_operator import DenseDiscreteBoundaryOperator

        product = 2 * self._operator_real
        self.assertIsInstance(product, DenseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_scalar_times_dense_is_correct(self):
        import numpy as np

        product = 2 * self._operator_real

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = 2 * self._operator_real * vec
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_product_of_dense_operators_is_dense_operator(self):
        from bempp.assembly.discrete_boundary_operator import DenseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator_real * self._operator_complex,
                              DenseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_dense_times_dense_is_correct(self):
        import numpy as np

        product = self._operator_real * self._operator_real

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = self._operator_real * (self._operator_real * vec)
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)


class TestSparseDiscreteBoundaryOperator(TestCase):
    """Test cases for sparse discrete operators."""

    def setUp(self):
        import bempp

        grid = bempp.grid_from_sphere(3)
        space = bempp.function_space(grid, "DP", 0)
        self._operator = bempp.operators.boundary.sparse.identity(space, space, space).weak_form()
        self._space = space

    def test_type(self):
        from bempp.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator, SparseDiscreteBoundaryOperator)

    def test_shape(self):
        expected = (self._space.global_dof_count, self._space.global_dof_count)
        actual = self._operator.shape
        self.assertEqual(expected, actual)

    def test_sum_of_sparse_operators_is_sparse(self):
        from bempp.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        sum_operator = self._operator + self._operator
        self.assertIsInstance(sum_operator, SparseDiscreteBoundaryOperator)

    def matrix_vector_product_of_sparse_operator_sum_is_correct(self):
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator + self._operator
        actual = sum_operator * vec
        expected = self._operator * vec + self._operator * vec

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_scalar_multiply_of_sparse_operator_is_sparse(self):
        from bempp.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        product = 2 * self._operator
        self.assertIsInstance(product, SparseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_scalar_times_sparse_is_correct(self):
        import numpy as np

        product = 2 * self._operator

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = 2 * self._operator * vec
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_product_of_sparse_operators_is_sparse_operator(self):
        from bempp.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator * self._operator,
                              SparseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_sparse_times_sparse_is_correct(self):
        import numpy as np

        product = self._operator * self._operator

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = self._operator * (self._operator * vec)
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)


class TestInverseSparseDiscreteBoundaryOperator(TestCase):
    def setUp(self):
        import bempp

        grid = bempp.shapes.sphere(h=0.2)
        self._space_const = bempp.function_space(grid, "DP", 0)
        self._space_lin = bempp.function_space(grid, "P", 1)

    def test_square_inverse(self):
        import bempp
        import numpy as np

        space = self._space_const
        op = bempp.operators.boundary.sparse.identity(space, space, space).weak_form()
        inverse_op = bempp.InverseSparseDiscreteBoundaryOperator(op)

        expected = np.eye(space.global_dof_count, space.global_dof_count)
        actual = op * inverse_op * expected

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    def test_inverse_of_thin_sparse_matrix(self):
        import bempp
        import numpy as np

        space_const = self._space_const
        space_lin = self._space_lin
        op = bempp.operators.boundary.sparse.identity(space_lin, space_const, space_const).weak_form()
        inverse_op = bempp.InverseSparseDiscreteBoundaryOperator(op)

        identity = np.eye(space_lin.global_dof_count, space_lin.global_dof_count)
        expected = np.eye(space_lin.global_dof_count, space_lin.global_dof_count)

        actual = inverse_op * (op * expected)

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    def test_inverse_of_thick_sparse_matrix(self):
        import bempp
        import numpy as np

        space_const = self._space_const
        space_lin = self._space_lin
        op = bempp.operators.boundary.sparse.identity(space_const, space_lin, space_lin).weak_form()
        inverse_op = bempp.InverseSparseDiscreteBoundaryOperator(op)

        identity = np.eye(space_lin.global_dof_count, space_lin.global_dof_count)
        expected = np.eye(space_lin.global_dof_count, space_lin.global_dof_count)

        actual = op * inverse_op * expected

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)


if __name__ == "__main__":
    from unittest import main

    main()
