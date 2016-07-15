from unittest import TestCase
import unittest
import bempp.api
import numpy as np


class TestHelmholtz(TestCase):

    def test_solve_low_order_helmholtz_dirichlet_problem_on_sphere(self):
        """Solve an exterior Helmholtz problem on the unit sphere."""

        grid = bempp.api.shapes.regular_sphere(4)
        lin_space = bempp.api.function_space(grid, "P", 1)
        wave_number = 1

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(wave_number * 1j * r) * (1j * wave_number * r - 1) / r**2

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(wave_number * 1j * r) / r

        grid_fun = bempp.api.GridFunction(lin_space, fun=fun)

        dlp = bempp.api.operators.boundary.helmholtz.double_layer(
            lin_space, lin_space, lin_space, wave_number)

        slp = bempp.api.operators.boundary.helmholtz.single_layer(
            lin_space, lin_space, lin_space, wave_number)

        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space, lin_space, lin_space)

        rhs = (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info = bempp.api.linalg.gmres(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(rel_error < 1E-2)

    def test_solve_high_order_helmholtz_dirichlet_problem_on_sphere(self):
        """Solve an exterior Helmholtz Dirichlet problem on the unit sphere using quadratic basis functions."""

        grid = bempp.api.shapes.regular_sphere(4)
        space = bempp.api.function_space(grid, "P", 2)

        parameters = bempp.api.common.global_parameters()
        parameters.quadrature.far.double_order = 3
        wave_number = 1

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(wave_number * 1j * r) * (1j * wave_number * r - 1) / r**2

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(wave_number * 1j * r) / r

        grid_fun = bempp.api.GridFunction(space, fun=fun)

        dlp = bempp.api.operators.boundary.helmholtz.double_layer(
            space, space, space, wave_number, parameters=parameters)

        slp = bempp.api.operators.boundary.helmholtz.single_layer(
            space, space, space, wave_number, parameters=parameters)

        ident = bempp.api.operators.boundary.sparse.identity(
            space, space, space, parameters=parameters)

        rhs = (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info = bempp.api.linalg.gmres(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(rel_error < 1E-2)

    def test_solve_helmholtz_neumann_problem_on_sphere(self):

        grid = bempp.api.shapes.regular_sphere(5)
        space = bempp.api.function_space(grid, "P", 2)

        parameters = bempp.api.common.global_parameters()
        parameters.quadrature.near.double_order = 6
        parameters.quadrature.medium.double_order = 3
        parameters.quadrature.far.double_order = 3
        wave_number = 1

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(wave_number * 1j * r) * \
                (1j * wave_number * r - 1) / r**2

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(wave_number * 1j * r) / r

        grid_fun = bempp.api.GridFunction(space, fun=fun)

        adlp = bempp.api.operators.boundary.helmholtz.adjoint_double_layer(
            space, space, space, wave_number, parameters=parameters)

        hyp = bempp.api.operators.boundary.helmholtz.hypersingular(
            space, space, space, wave_number, parameters=parameters)

        ident = bempp.api.operators.boundary.sparse.identity(
            space, space, space, parameters=parameters)

        rhs = (-.5 * ident - adlp) * grid_fun

        sol_grid_fun, info = bempp.api.linalg.gmres(hyp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(rel_error < 1E-2)

    def test_operator_preconditioned_dirichlet_solve_dual(self):
        """Solve an operator preconditioned Helmholtz Dirichlet problem with dual spaces."""

        grid = bempp.api.shapes.regular_sphere(3)

        k = 1

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(1j * k * r) / r**2 * (1j * k * r - 1)

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(1j * k * r)

        slp, hyp = bempp.api.operators.boundary.helmholtz.single_layer_and_hypersingular_pair(
            grid, k, spaces='dual')

        const_space = slp.domain
        lin_space = slp.range

        grid_fun = bempp.api.GridFunction(
            lin_space, dual_space=const_space, fun=fun)

        dlp = bempp.api.operators.boundary.helmholtz.double_layer(
            lin_space, lin_space, const_space, k)

        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space, lin_space, const_space)

        rhs = hyp * (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info, res = bempp.api.linalg.gmres(
            hyp * slp, rhs, return_residuals=True)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(
            rel_error < 2E-2, msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))
        self.assertTrue(len(
            res) < 9, msg="Needed {0} iterations to solve system. Expected not more than 8 iterations.".format(len(res)))

    def test_operator_preconditioned_dirichlet_solve_lin(self):
        """Solve an operator preconditioned Helmholtz Dirichlet problem with linear spaces."""

        grid = bempp.api.shapes.regular_sphere(3)

        k = 1

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(1j * k * r) / r**2 * (1j * k * r - 1)

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(1j * k * r)

        slp, hyp = bempp.api.operators.boundary.helmholtz.single_layer_and_hypersingular_pair(
            grid, k, spaces='linear')

        const_space = slp.domain
        lin_space = slp.range

        grid_fun = bempp.api.GridFunction(
            lin_space, dual_space=const_space, fun=fun)

        dlp = bempp.api.operators.boundary.helmholtz.double_layer(
            lin_space, lin_space, const_space, k)

        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space, lin_space, const_space)

        rhs = hyp * (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info, res = bempp.api.linalg.gmres(
            hyp * slp, rhs, return_residuals=True)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(
            rel_error < 2E-2, msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))
        self.assertTrue(len(
            res) < 10, msg="Needed {0} iterations to solve system. Expected not more than 9 iterations.".format(len(res)))

    def test_helmholtz_potentials(self):
        """Test the Helmholtz potential operators."""

        grid = bempp.api.shapes.regular_sphere(4)

        space = bempp.api.function_space(grid, "P", 1)

        def dirichlet_fun(x):
            return np.exp(1j * x[0])

        def dirichlet_data(x, n, domain_index, res):
            res[0] = dirichlet_fun(x)

        def neumann_data(x, n, domain_index, res):
            res[0] = 1j * np.exp(1j * x[0]) * n[0]

        dirichlet_grid_fun = bempp.api.GridFunction(space, fun=dirichlet_data)
        neumann_grid_fun = bempp.api.GridFunction(space, fun=neumann_data)

        point = np.array([[0, 0.2, 0.3]]).T

        sl = bempp.api.operators.potential.helmholtz.single_layer(
            space, point, 1)
        dl = bempp.api.operators.potential.helmholtz.double_layer(
            space, point, 1)

        actual = (sl * neumann_grid_fun - dl * dirichlet_grid_fun)[0, 0]
        expected = dirichlet_fun(point[:, 0])

        rel_error = np.abs(actual - expected) / np.abs(expected)
        self.assertTrue(
            rel_error < 1E-6, msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))

if __name__ == "__main__":

    from unittest import main
    main()
