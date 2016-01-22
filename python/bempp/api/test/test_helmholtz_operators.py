from unittest import TestCase
import unittest
import bempp.api
import numpy as np


class TestHelmholtzOperators(TestCase):

    def test_solve_low_order_helmholtz_dirichlet_problem_on_sphere(self):
        """Solve an exterior Laplace problem on the unit sphere."""

        grid = bempp.api.shapes.sphere(h=0.1)
        const_space = bempp.api.function_space(grid, "DP", 0)
        lin_space = bempp.api.function_space(grid, "P", 1)
        wave_number = 1

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(wave_number * 1j * r) * ( 1j * wave_number * r - 1) / r**2

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(wave_number * 1j * r) / r

        grid_fun = bempp.api.GridFunction(lin_space, dual_space = const_space, fun=fun)

        dlp = bempp.api.operators.boundary.helmholtz.double_layer(
                lin_space, lin_space, const_space, wave_number)

        slp = bempp.api.operators.boundary.helmholtz.single_layer(
                const_space, lin_space, const_space, wave_number)

        ident = bempp.api.operators.boundary.sparse.identity(
                lin_space, lin_space, const_space)

        rhs = (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info = bempp.api.linalg.gmres(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(rel_error < 1E-2)


    def test_solve_high_order_helmholtz_dirichlet_problem_on_sphere(self):
        """Solve an exterior Laplace Dirichlet problem on the unit sphere using quadratic basis functions."""

        grid = bempp.api.shapes.sphere(h=0.1)
        space = bempp.api.function_space(grid, "P", 2)

        parameters = bempp.api.common.global_parameters()
        parameters.quadrature.far.double_order = 3
        wave_number = 1

        def sol_fun(x):
            r = np.linalg.norm(x)
            return np.exp(wave_number * 1j * r) * ( 1j * wave_number * r - 1) / r**2

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

        grid = bempp.api.shapes.sphere(h=0.1)
        space = bempp.api.function_space(grid, "P", 2)

        parameters = bempp.api.common.global_parameters()
        parameters.quadrature.far.double_order = 3
        wave_number = 1

        def fun(x, n, domain_index, res):
            r = np.linalg.norm(x)
            res[0] = np.exp(wave_number * 1j * r) * ( 1j * wave_number * r - 1) / r**2

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


if __name__ == "__main__":

    from unittest import main
    main()


