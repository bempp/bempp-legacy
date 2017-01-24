"""Test discrete boundary operators."""

from unittest import TestCase
import unittest
import bempp.api

#pylint: disable=invalid-name

class TestGeneralNonlocalDiscreteBoundaryOperator(TestCase):
    """Test cases for general discrete operators."""

    def setUp(self):
        """Setup the test cases."""

        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'hmat'
        self._parameters = parameters
        grid = bempp.api.shapes.regular_sphere(3)
        space = bempp.api.function_space(grid, "DP", 0)
        self._operator_real = \
            bempp.api.operators.boundary.laplace.single_layer(
                space, space, space, parameters=parameters,
                use_projection_spaces=False).weak_form()
        self._operator_complex = \
            bempp.api.operators.boundary.helmholtz.single_layer(
                space, space, space, 1,
                parameters=parameters,
                use_projection_spaces=False).weak_form()
        self._space = space

    def test_type(self):
        """Type of discrete operator is correct."""
        from bempp.api.assembly.discrete_boundary_operator import \
            GeneralNonlocalDiscreteBoundaryOperator

        self.assertIsInstance(self._operator_real,
                              GeneralNonlocalDiscreteBoundaryOperator)

    def test_shape(self):
        """Discrete operator has correct shape."""
        expected = (self._space.global_dof_count, self._space.global_dof_count)
        actual = self._operator_real.shape
        self.assertEqual(expected, actual)

    def test_sum_of_real_and_complex_has_complex_type(self):
        """Sum of real and complex operator has complex type."""
        sum_operator = self._operator_real + self._operator_complex
        self.assertEqual(sum_operator.dtype, 'complex128')

    def matrix_vector_product_of_operator_sum_is_correct(self):
        """Matrix vector product of sum is correct."""
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        sum_operator = self._operator_real + self._operator_complex
        actual = sum_operator * vec
        expected = self._operator_real * vec + self._operator_complex * vec
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)


class TestDenseDiscreteBoundaryOperator(TestCase):
    """Test cases for dense discrete operators."""

    def setUp(self):
        """Setup the test cases."""
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        self._parameters = parameters
        grid = bempp.api.shapes.regular_sphere(3)
        space = bempp.api.function_space(grid, "DP", 0)
        self._operator_real = \
            bempp.api.operators.boundary.laplace.single_layer(
                space, space, space, parameters=parameters,
                use_projection_spaces=False).weak_form()
        self._operator_complex = \
            bempp.api.operators.boundary.helmholtz.single_layer(
                space, space, space, 1, parameters=parameters,
                use_projection_spaces=False).weak_form()
        self._space = space

    def test_type(self):
        """Type is correct."""
        from bempp.api.assembly.discrete_boundary_operator import \
            DenseDiscreteBoundaryOperator

        self.assertIsInstance(
            self._operator_real, DenseDiscreteBoundaryOperator)

    def test_shape(self):
        """Shape is correct."""
        expected = (self._space.global_dof_count, self._space.global_dof_count)
        actual = self._operator_real.shape
        self.assertEqual(expected, actual)

    def test_sum_of_real_and_complex_has_complex_type(self):
        """Sum of real and complex dense discrete operator has complex type."""
        sum_operator = self._operator_real + self._operator_complex
        self.assertEqual(sum_operator.dtype, 'complex128')

    def test_sum_of_dense_operators_is_dense(self):
        """Some of dense operators is dense."""
        from bempp.api.assembly.discrete_boundary_operator import \
            DenseDiscreteBoundaryOperator

        sum_operator = self._operator_real + self._operator_complex
        self.assertIsInstance(sum_operator, DenseDiscreteBoundaryOperator)

    def matrix_vector_product_of_operator_sum_is_correct(self):
        """Matrix vector product of sum is correct for dense discrete op."""
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator_real + self._operator_complex
        actual = sum_operator * vec
        expected = self._operator_real * vec + self._operator_complex * vec

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_scalar_multiply_of_dense_operator_is_dense(self):
        """Scalar Multiply of dense operator is dense."""
        from bempp.api.assembly.discrete_boundary_operator import \
            DenseDiscreteBoundaryOperator
        product = 2 * self._operator_real
        self.assertIsInstance(product, DenseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_scalar_times_dense_is_correct(self):
        """Matvec of scalar times dense operator is correct."""
        import numpy as np

        product = 2 * self._operator_real
        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = 2 * self._operator_real * vec
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_product_of_dense_operators_is_dense_operator(self):
        """Product of dense operators is dense."""
        from bempp.api.assembly.discrete_boundary_operator import \
            DenseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator_real * self._operator_complex,
                              DenseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_dense_times_dense_is_correct(self):
        """Matvec of dense times dense is correct."""
        import numpy as np

        product = self._operator_real * self._operator_real

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = self._operator_real * (self._operator_real * vec)
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)


class TestSparseDiscreteBoundaryOperator(TestCase):
    """Test cases for sparse discrete operators."""

    def setUp(self):
        """Setup the test cases."""
        grid = bempp.api.shapes.regular_sphere(3)
        space = bempp.api.function_space(grid, "DP", 0)
        self._operator = bempp.api.operators.boundary.sparse.identity(
            space, space, space).weak_form()
        self._space = space

    def test_type(self):
        """Type is correct."""
        from bempp.api.assembly.discrete_boundary_operator import \
            SparseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator, SparseDiscreteBoundaryOperator)

    def test_shape(self):
        """Shape is correct."""
        expected = (self._space.global_dof_count, self._space.global_dof_count)
        actual = self._operator.shape
        self.assertEqual(expected, actual)

    def test_sum_of_sparse_operators_is_sparse(self):
        """Sum of sparse discrete operators is sparse."""
        from bempp.api.assembly.discrete_boundary_operator import \
            SparseDiscreteBoundaryOperator

        sum_operator = self._operator + self._operator
        self.assertIsInstance(sum_operator, SparseDiscreteBoundaryOperator)

    def matrix_vector_product_of_sparse_operator_sum_is_correct(self):
        """Matvec of sparse operator sum is correct."""
        import numpy as np

        vec = np.ones(self._space.global_dof_count, dtype='float64')

        sum_operator = self._operator + self._operator
        actual = sum_operator * vec
        expected = self._operator * vec + self._operator * vec

        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_scalar_multiply_of_sparse_operator_is_sparse(self):
        """Scalar multiply of sparse operator is sparse."""
        from bempp.api.assembly.discrete_boundary_operator import \
            SparseDiscreteBoundaryOperator

        product = 2 * self._operator
        self.assertIsInstance(product, SparseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_scalar_times_sparse_is_correct(self):
        """Matvec of scalar times sparse is correct."""
        import numpy as np

        product = 2 * self._operator

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = 2 * self._operator * vec
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)

    def test_product_of_sparse_operators_is_sparse_operator(self):
        """Product of sparse operators is sparse operator."""
        from bempp.api.assembly.discrete_boundary_operator import \
            SparseDiscreteBoundaryOperator

        self.assertIsInstance(self._operator * self._operator,
                              SparseDiscreteBoundaryOperator)

    def test_matrix_vector_product_of_sparse_times_sparse_is_correct(self):
        """Matvec of sparse operator is correct."""
        import numpy as np

        product = self._operator * self._operator

        vec = np.ones(self._space.global_dof_count, dtype='float64')
        actual = product * vec
        expected = self._operator * (self._operator * vec)
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)


class TestInverseSparseDiscreteBoundaryOperator(TestCase):
    """Test cases for inverse sparse discrete operators."""

    requiresgmsh = unittest.skipIf(
        bempp.api.GMSH_PATH is None, reason="Needs GMSH")

    def setUp(self):
        """Setup the test cases."""
        grid = bempp.api.shapes.sphere(h=0.2)
        self._space_const = bempp.api.function_space(grid, "DP", 0)
        self._space_lin = bempp.api.function_space(grid, "P", 1)

    @requiresgmsh
    def test_square_inverse(self):
        """Test the inverse of a square sparse matrix."""
        import numpy as np

        space = self._space_const
        op = bempp.api.operators.boundary.sparse.identity(
            space, space, space).weak_form()
        inverse_op = bempp.api.InverseSparseDiscreteBoundaryOperator(op)

        expected = np.eye(space.global_dof_count, space.global_dof_count)
        actual = op * inverse_op * expected

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    @requiresgmsh
    def test_inverse_of_thin_sparse_matrix(self):
        """Test the pseudo-inverse of a thin sparse matrix."""
        import numpy as np

        space_const = self._space_const
        space_lin = self._space_lin
        op = bempp.api.operators.boundary.sparse.identity(
            space_lin, space_const, space_const).weak_form()
        inverse_op = bempp.api.InverseSparseDiscreteBoundaryOperator(op)

        expected = np.eye(space_lin.global_dof_count,
                          space_lin.global_dof_count)

        actual = inverse_op * (op * expected)

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    @requiresgmsh
    def test_inverse_of_thick_sparse_matrix(self):
        """Pseudo-Inverse of thick sparse matrix."""
        import numpy as np

        space_const = self._space_const
        space_lin = self._space_lin
        op = bempp.api.operators.boundary.sparse.identity(
            space_const, space_lin, space_lin).weak_form()
        inverse_op = bempp.api.InverseSparseDiscreteBoundaryOperator(op)

        expected = np.eye(space_lin.global_dof_count,
                          space_lin.global_dof_count)

        actual = op * inverse_op * expected

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)


class TestZeroDiscreteBoundaryOperator(TestCase):
    """Test cases for zero discrete operator."""

    def setUp(self):
        """Setup test cases."""
        from bempp.api.assembly.discrete_boundary_operator import \
            ZeroDiscreteBoundaryOperator

        self._M = 10
        self._N = 5
        self._op = ZeroDiscreteBoundaryOperator(self._M, self._N)

    def test_shape(self):
        """Shape is correct."""

        self.assertEqual(self._op.shape, (self._M, self._N))

    def test_multiply_with_vector(self):
        """Multiplication with vector is correct."""
        import numpy as np
        x = np.ones((self._N, 1), dtype='float64')
        res = self._op * x
        self.assertEqual(res.shape, (self._M, 1),
                         "Multiply with array with ndim = 2.")

        res = self._op * x.squeeze()
        self.assertEqual(res.shape, (self._M, ),
                         "Multiply with array with ndim = 1.")


if __name__ == "__main__":
    from unittest import main
    main()
