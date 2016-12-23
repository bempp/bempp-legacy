"""Test cases for Grid Functions."""

from unittest import TestCase
import bempp.api

#pylint: disable=invalid-name

class TestGridFunction(TestCase):
    """Test class for the GridFunction object."""

    def setUp(self):
        """Setup test cases."""
        grid = bempp.api.shapes.regular_sphere(3)
        self._space = bempp.api.function_space(grid, "P", 1)

    def test_initialize_from_coefficients(self):
        """Initialize grid function from coefficients."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(
            self._space, coefficients=coefficients)

        self.assertAlmostEqual(np.linalg.norm(
            coefficients - grid_fun.coefficients), 0)

    def test_set_coefficients(self):
        """Set coefficients of grid function."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(
            self._space, coefficients=coefficients)
        grid_fun.coefficients *= 2
        self.assertAlmostEqual(np.linalg.norm(
            2 * coefficients - grid_fun.coefficients), 0)

    def test_initialize_from_real_function(self):
        """Initialize from real function."""
        import numpy as np

        #pylint: disable=unused-argument
        def fun(point, normal, domain_index, res):
            """Discretization function."""
            res[0] = point[0]

        grid_fun = bempp.api.GridFunction(self._space, fun=fun)
        actual = grid_fun.coefficients
        expected = self._space.global_dof_interpolation_points[0, :]
        self.assertAlmostEqual(np.linalg.norm(actual - expected), 0)
        self.assertEqual(grid_fun.dtype, 'float64')

    def test_initialize_from_complex_function(self):
        """Initialize from complex function."""
        import numpy as np

        #pylint: disable=unused-argument
        def fun(point, normal, domain_index, res):
            """Discretization function."""
            res[0] = 1j * point[0]

        grid_fun = bempp.api.GridFunction(self._space, fun=fun)
        actual = grid_fun.coefficients
        expected = self._space.global_dof_interpolation_points[0, :]
        self.assertAlmostEqual(np.linalg.norm(actual - 1j * expected), 0)
        self.assertEqual(grid_fun.dtype, 'complex128')

    def test_initialize_from_projections(self):
        """Initialize from projections."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(
            self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        self.assertAlmostEqual(np.linalg.norm(
            coefficients - grid_fun.coefficients), 0)

    def test_add_two_grid_functions(self):
        """Add two grid functions."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(
            self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = grid_fun + grid_fun

        expected = 2 * grid_fun.coefficients
        actual = grid_fun2.coefficients

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    def test_scalar_multiply_grid_function(self):
        """Multiply grid function with a scalar."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(
            self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = 2 * grid_fun

        expected = 2 * grid_fun.coefficients
        actual = grid_fun2.coefficients

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    def test_scalar_divide_grid_function(self):
        """Divide grid function by a scalar."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(
            self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = grid_fun / 2

        expected = grid_fun.coefficients / 2
        actual = grid_fun2.coefficients

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    def test_negate_grid_function(self):
        """Negate grid function."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        identity = bempp.api.operators.boundary.sparse.identity(
            self._space, self._space, self._space).weak_form()
        projections = identity * coefficients

        grid_fun = bempp.api.GridFunction(self._space, projections=projections)

        grid_fun2 = -grid_fun

        expected = -grid_fun.coefficients
        actual = grid_fun2.coefficients

        self.assertAlmostEqual(np.linalg.norm(expected - actual), 0)

    def test_l2_norm(self):
        """L2 norm of a grid function."""
        import numpy as np

        space = self._space
        n = space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(space, coefficients=coefficients)

        ident = bempp.api.operators.boundary.sparse.identity(
            space, space, space)

        expected = np.sqrt(
            np.dot(coefficients, ident.weak_form() * coefficients))
        actual = grid_fun.l2_norm()

        self.assertAlmostEqual(expected, actual, 9)

    def test_sum_of_l2_norms_on_elements_is_equal_to_full_l2_norm(self):
        """Sum of L2 norms equals full L2 norm."""
        import numpy as np

        n = self._space.global_dof_count
        coefficients = np.ones(n)
        grid_fun = bempp.api.GridFunction(
            self._space, coefficients=coefficients)

        l2_sum = 0

        for element in grid_fun.space.grid.leaf_view.entity_iterator(0):
            l2_sum += grid_fun.l2_norm(element)**2

        self.assertAlmostEqual(np.sqrt(l2_sum), grid_fun.l2_norm())

    #pylint: disable=too-many-locals
    def test_surface_grad_is_correct(self):
        """Compute the correct surface gradient."""

        import numpy as np
        parameters = bempp.api.common.global_parameters()
        parameters.quadrature.far.single_order = 8
        grid = bempp.api.shapes.regular_sphere(3)
        space = bempp.api.function_space(grid, "P", 4)

        elements = list(grid.leaf_view.entity_iterator(0))

        #pylint: disable=unused-argument
        def fun(point, normal, domain_index, result):
            """Discretization function."""
            result[0] = np.sum(point**2)

        grid_fun = bempp.api.GridFunction(
            space, fun=fun, parameters=parameters)

        element = elements[0]
        geometry = element.geometry
        local_coordinate = np.array([[0], [0]])
        global_coordinate = geometry.local2global(local_coordinate)
        normal = geometry.normals(local_coordinate)
        grad = 2 * global_coordinate[:, 0]

        expected = grad - normal[:, 0] * np.dot(grad, normal[:, 0])
        actual = grid_fun.evaluate_surface_gradient(
            element, local_coordinate)[:, 0]
        rel_diff = np.linalg.norm(expected - actual) / np.linalg.norm(actual)

        self.assertAlmostEqual(rel_diff, 0)


if __name__ == "__main__":
    from unittest import main

    main()
