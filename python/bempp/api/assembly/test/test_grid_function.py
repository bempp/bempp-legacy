"""Test cases for Grid Functions."""

from unittest import TestCase
import bempp.api


class TestGridFunction(TestCase):
    """Test class for the GridFunction object."""

    def setUp(self):
        grid = bempp.api.shapes.regular_sphere(3)
        self._space = bempp.api.function_space(grid, "P", 1)

    def test_initialize_from_coefficients(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(self._space, coefficients=coefficients)

        self.assertAlmostEquals(np.linalg.norm(coefficients - grid_fun.coefficients), 0)

    def test_set_coefficients(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(self._space, coefficients=coefficients)
        grid_fun.coefficients *= 2
        self.assertAlmostEquals(np.linalg.norm(2 * coefficients - grid_fun.coefficients), 0)

    def test_initialize_from_real_function(self):
        import numpy as np

        def fun(x, n, d, res):
            res[0] = x[0]

        grid_fun = bempp.api.GridFunction(self._space, fun=fun)
        actual = grid_fun.coefficients
        expected = self._space.global_dof_interpolation_points[0, :]
        self.assertAlmostEquals(np.linalg.norm(actual - expected), 0)
        self.assertEqual(grid_fun.dtype, 'float64')

    def test_initialize_from_complex_function(self):
        import numpy as np

        def fun(x, n, d, res):
            res[0] = 1j * x[0]

        grid_fun = bempp.api.GridFunction(self._space, fun=fun)
        actual = grid_fun.coefficients
        expected = self._space.global_dof_interpolation_points[0, :]
        self.assertAlmostEquals(np.linalg.norm(actual - 1j * expected), 0)
        self.assertEqual(grid_fun.dtype, 'complex128')

    def test_initialize_from_projections(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        self.assertAlmostEquals(np.linalg.norm(coefficients - grid_fun.coefficients), 0)

    def test_add_two_grid_functions(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = grid_fun + grid_fun

        expected = 2 * grid_fun.coefficients
        actual = grid_fun2.coefficients

        self.assertAlmostEquals(np.linalg.norm(expected - actual), 0)

    def test_scalar_multiply_grid_function(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = 2 * grid_fun

        expected = 2 * grid_fun.coefficients
        actual = grid_fun2.coefficients

        self.assertAlmostEquals(np.linalg.norm(expected - actual), 0)

    def test_scalar_divide_grid_function(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = grid_fun/2

        expected = grid_fun.coefficients/2
        actual = grid_fun2.coefficients

        self.assertAlmostEquals(np.linalg.norm(expected - actual), 0)

    def test_negate_grid_function(self):
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = -grid_fun

        expected = -grid_fun.coefficients
        actual = grid_fun2.coefficients

        self.assertAlmostEquals(np.linalg.norm(expected - actual), 0)

    def test_sum_of_l2_norms_on_elements_is_equal_to_full_l2_norm(self):

        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(self._space, coefficients=coefficients)

        sum = 0

        for element in grid_fun.space.grid.leaf_view.entity_iterator(0):
            sum += grid_fun.l2_norm(element)**2

        self.assertAlmostEqual(np.sqrt(sum), grid_fun.l2_norm())


if __name__ == "__main__":
    from unittest import main

    main()
