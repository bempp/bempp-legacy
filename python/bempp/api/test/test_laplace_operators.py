from unittest import TestCase
import bempp.api
import numpy as np


class TestLaplaceOperators(TestCase):

    def test_solve_low_order_laplace_dirichlet_problem_on_sphere(self):
        """Solve an exterior Laplace problem on the unit sphere."""

        grid = bempp.api.shapes.sphere(h=0.1)
        const_space = bempp.api.function_space(grid, "DP", 0)
        lin_space = bempp.api.function_space(grid, "P", 1)

        def sol_fun(x):
            return -1./np.linalg.norm(x)**2

        def fun(x, n, domain_index, res):
            res[0] = 1

        grid_fun = bempp.api.GridFunction(lin_space, dual_space = const_space, fun=fun)

        dlp = bempp.api.operators.boundary.laplace.double_layer(
                lin_space, lin_space, const_space)

        slp = bempp.api.operators.boundary.laplace.single_layer(
                const_space, lin_space, const_space)

        ident = bempp.api.operators.boundary.sparse.identity(
                lin_space, lin_space, const_space)

        rhs = (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info = bempp.api.linalg.cg(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(rel_error < 1E-2)


    def test_solve_high_order_laplace_dirichlet_problem_on_sphere(self):
        """Solve an exterior Laplace Dirichlet problem on the unit sphere using quadratic basis functions."""

        grid = bempp.api.shapes.sphere(h=0.1)
        space = bempp.api.function_space(grid, "P", 2)

        parameters = bempp.api.common.global_parameters()
        parameters.quadrature.far.double_order = 3

        def sol_fun(x):
            return -1./np.linalg.norm(x)**2

        def fun(x, n, domain_index, res):
            res[0] = 1

        grid_fun = bempp.api.GridFunction(space, fun=fun)

        dlp = bempp.api.operators.boundary.laplace.double_layer(
                space, space, space, parameters=parameters)

        slp = bempp.api.operators.boundary.laplace.single_layer(
                space, space, space, parameters=parameters)

        ident = bempp.api.operators.boundary.sparse.identity(
                space, space, space, parameters=parameters)

        rhs = (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, info = bempp.api.linalg.cg(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(rel_error < 1E-2)

if __name__ == "__main__":

    from unittest import main
    main()


