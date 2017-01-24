"""Various validation tests for Laplace problems."""
from unittest import TestCase
import bempp.api
import numpy as np

#pylint: disable=invalid-name
#pylint: disable=unused-argument
#pylint: disable=missing-docstring
#pylint: disable=too-many-locals


class TestLaplace(TestCase):
    """Laplace test cases."""

    def test_solve_low_order_laplace_dirichlet_problem_on_sphere(self):
        """Solve an exterior Laplace problem on the unit sphere."""

        grid = bempp.api.shapes.regular_sphere(4)
        lin_space = bempp.api.function_space(grid, "P", 1)

        def sol_fun(x):
            return -1. / np.linalg.norm(x)**2

        def fun(x, n, domain_index, res):
            res[0] = 1

        grid_fun = bempp.api.GridFunction(lin_space, fun=fun)

        dlp = bempp.api.operators.boundary.laplace.double_layer(
            lin_space, lin_space, lin_space)

        slp = bempp.api.operators.boundary.laplace.single_layer(
            lin_space, lin_space, lin_space)

        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space, lin_space, lin_space)

        rhs = (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, _ = bempp.api.linalg.cg(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(
            rel_error < 1E-2,
            msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))

    def test_solve_high_order_laplace_dirichlet_problem_on_sphere(self):
        """Exterior Dirichlet problem with quadratic basis functions."""

        grid = bempp.api.shapes.regular_sphere(4)
        space = bempp.api.function_space(grid, "P", 2)

        parameters = bempp.api.common.global_parameters()
        parameters.hmat.eps = 1E-5
        parameters.quadrature.medium.double_order = 4
        parameters.quadrature.far.double_order = 4

        def sol_fun(x):
            return -1. / np.linalg.norm(x)**2

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

        sol_grid_fun, _ = bempp.api.linalg.cg(slp, rhs)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(
            rel_error < 1E-2,
            msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))

    def test_null_space_of_hypersingular_operator(self):
        """Test the null space of a hypersingular operator."""

        grid = bempp.api.shapes.regular_sphere(4)
        lin_space = bempp.api.function_space(grid, "P", 1)

        def fun(x, n, domain_index, res):
            res[0] = 1

        grid_fun = bempp.api.GridFunction(lin_space, fun=fun)

        hyp = bempp.api.operators.boundary.laplace.hypersingular(
            lin_space, lin_space, lin_space)

        res = hyp * grid_fun

        rel_error = res.l2_norm() / grid_fun.l2_norm()
        self.assertTrue(
            rel_error < 3E-3,
            msg="Actual error: {0}. Expected error: 3E-3".format(rel_error))

    def test_mixed_laplace_problem(self):
        """Mixed Laplace problem with higher order basis functions."""

        grid = bempp.api.shapes.cube()
        dirichlet_segments = [1, 3]
        neumann_segments = [2, 4, 5, 6]

        parameters = bempp.api.common.global_parameters()

        parameters.quadrature.medium.double_order = 4
        parameters.quadrature.far.double_order = 4

        order_neumann = 1
        order_dirichlet = 2

        global_dirichlet_space = bempp.api.function_space(
            grid, "P", order_dirichlet)

        neumann_space_dirichlet_segment = bempp.api.function_space(
            grid, "DP", order_neumann,
            domains=dirichlet_segments, closed=True,
            element_on_segment=True)
        neumann_space_neumann_segment = bempp.api.function_space(
            grid, "DP", order_neumann,
            domains=neumann_segments, closed=False,
            element_on_segment=True,
            reference_point_on_segment=False)
        dirichlet_space_dirichlet_segment = bempp.api.function_space(
            grid, "P", order_dirichlet,
            domains=dirichlet_segments, closed=True)
        dirichlet_space_neumann_segment = bempp.api.function_space(
            grid, "P", order_dirichlet,
            domains=neumann_segments, closed=False)
        dual_dirichlet_space = bempp.api.function_space(
            grid, "P", order_dirichlet, domains=dirichlet_segments,
            closed=True, strictly_on_segment=True)
        slp_DD = bempp.api.operators.boundary.laplace.single_layer(
            neumann_space_dirichlet_segment,
            dirichlet_space_dirichlet_segment,
            neumann_space_dirichlet_segment,
            parameters=parameters)
        dlp_DN = bempp.api.operators.boundary.laplace.double_layer(
            dirichlet_space_neumann_segment,
            dirichlet_space_dirichlet_segment,
            neumann_space_dirichlet_segment,
            parameters=parameters)
        adlp_ND = bempp.api.operators.boundary.laplace.adjoint_double_layer(
            neumann_space_dirichlet_segment,
            neumann_space_neumann_segment,
            dirichlet_space_neumann_segment,
            parameters=parameters)
        hyp_NN = bempp.api.operators.boundary.laplace.hypersingular(
            dirichlet_space_neumann_segment,
            neumann_space_neumann_segment,
            dirichlet_space_neumann_segment,
            parameters=parameters)
        slp_DN = bempp.api.operators.boundary.laplace.single_layer(
            neumann_space_neumann_segment,
            dirichlet_space_dirichlet_segment,
            neumann_space_dirichlet_segment,
            parameters=parameters)
        dlp_DD = bempp.api.operators.boundary.laplace.double_layer(
            dirichlet_space_dirichlet_segment,
            dirichlet_space_dirichlet_segment,
            neumann_space_dirichlet_segment,
            parameters=parameters)
        id_DD = bempp.api.operators.boundary.sparse.identity(
            dirichlet_space_dirichlet_segment,
            dirichlet_space_dirichlet_segment,
            neumann_space_dirichlet_segment,
            parameters=parameters)
        adlp_NN = bempp.api.operators.boundary.laplace.adjoint_double_layer(
            neumann_space_neumann_segment,
            neumann_space_neumann_segment,
            dirichlet_space_neumann_segment,
            parameters=parameters)
        id_NN = bempp.api.operators.boundary.sparse.identity(
            neumann_space_neumann_segment,
            neumann_space_neumann_segment,
            dirichlet_space_neumann_segment,
            parameters=parameters)
        hyp_ND = bempp.api.operators.boundary.laplace.hypersingular(
            dirichlet_space_dirichlet_segment,
            neumann_space_neumann_segment,
            dirichlet_space_neumann_segment,
            parameters=parameters)
        blocked = bempp.api.BlockedOperator(2, 2)
        blocked[0, 0] = slp_DD
        blocked[0, 1] = -dlp_DN
        blocked[1, 0] = adlp_ND
        blocked[1, 1] = hyp_NN

        def dirichlet_data_fun(x):
            return np.exp(x[0]) * np.sin(x[1])

        def dirichlet_data(x, n, domain_index, res):
            res[0] = dirichlet_data_fun(x)

        def neumann_data(x, n, domain_index, res):
            grad = np.array([np.exp(x[0]) * np.sin(x[1]),
                             np.exp(x[0]) * np.cos(x[1]), 0])
            res[0] = np.dot(grad, n)

        dirichlet_grid_fun = bempp.api.GridFunction(
            dirichlet_space_dirichlet_segment,
            fun=dirichlet_data,
            dual_space=dual_dirichlet_space)

        neumann_grid_fun = bempp.api.GridFunction(
            neumann_space_neumann_segment,
            fun=neumann_data,
            dual_space=dirichlet_space_neumann_segment)

        rhs_fun1 = (.5 * id_DD + dlp_DD) * dirichlet_grid_fun - \
            slp_DN * neumann_grid_fun
        rhs_fun2 = -hyp_ND * dirichlet_grid_fun + \
            (.5 * id_NN - adlp_NN) * neumann_grid_fun

        lhs = blocked.weak_form()
        rhs = np.hstack([rhs_fun1.projections(neumann_space_dirichlet_segment),
                         rhs_fun2.projections(dirichlet_space_neumann_segment)])

        from scipy.sparse.linalg import gmres
        x, _ = gmres(lhs, rhs)

        nx0 = neumann_space_dirichlet_segment.global_dof_count
        dirichlet_solution = bempp.api.GridFunction(
            dirichlet_space_neumann_segment, coefficients=x[nx0:])

        dirichlet_imbedding_dirichlet_segment = \
            bempp.api.operators.boundary.sparse.identity(
                dirichlet_space_dirichlet_segment,
                global_dirichlet_space,
                global_dirichlet_space)

        dirichlet_imbedding_neumann_segment = \
            bempp.api.operators.boundary.sparse.identity(
                dirichlet_space_neumann_segment,
                global_dirichlet_space,
                global_dirichlet_space)

        dirichlet = (dirichlet_imbedding_dirichlet_segment *
                     dirichlet_grid_fun +
                     dirichlet_imbedding_neumann_segment * dirichlet_solution)

        rel_error = dirichlet.relative_error(dirichlet_data_fun)
        self.assertTrue(
            rel_error < 2E-3,
            msg="Actual error: {0}. Expected error: 2E-3".format(rel_error))

    def test_operator_preconditioned_dirichlet_solve_dual(self):
        """Operator preconditioned Dirichlet problem with dual spaces."""

        grid = bempp.api.shapes.regular_sphere(4)

        def sol_fun(x):
            return -1. / np.linalg.norm(x)**2

        def fun(x, n, domain_index, res):
            res[0] = 1

        from bempp.api.operators.boundary.laplace import \
            single_layer_and_hypersingular_pair

        slp, hyp = single_layer_and_hypersingular_pair(
            grid, spaces='dual', stabilization_factor=1)

        const_space = slp.domain
        lin_space = slp.range

        grid_fun = bempp.api.GridFunction(
            lin_space, dual_space=const_space, fun=fun)

        dlp = bempp.api.operators.boundary.laplace.double_layer(
            lin_space, lin_space, const_space)

        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space, lin_space, const_space)

        rhs = hyp * (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, _, res = bempp.api.linalg.cg(
            hyp * slp, rhs, return_residuals=True)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(
            rel_error < 1E-2,
            msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))
        self.assertTrue(
            len(res) < 9,
            msg="Needed {0} iterations to solve system. ".format(len(res)) +
            "Expected not more than 8 iterations.")

    def test_operator_preconditioned_dirichlet_solve_lin(self):
        """Operator preconditioned Dirichlet problem with linear spaces."""

        grid = bempp.api.shapes.regular_sphere(4)

        def sol_fun(x):
            return -1. / np.linalg.norm(x)**2

        def fun(x, n, domain_index, res):
            res[0] = 1

        from bempp.api.operators.boundary.laplace import \
            single_layer_and_hypersingular_pair

        slp, hyp = single_layer_and_hypersingular_pair(
            grid, spaces='linear', stabilization_factor=1)

        const_space = slp.domain
        lin_space = slp.range

        grid_fun = bempp.api.GridFunction(
            lin_space, dual_space=const_space, fun=fun)

        dlp = bempp.api.operators.boundary.laplace.double_layer(
            lin_space, lin_space, const_space)

        ident = bempp.api.operators.boundary.sparse.identity(
            lin_space, lin_space, const_space)

        rhs = hyp * (-.5 * ident + dlp) * grid_fun

        sol_grid_fun, _, res = bempp.api.linalg.cg(
            hyp * slp, rhs, return_residuals=True)

        rel_error = sol_grid_fun.relative_error(sol_fun)
        self.assertTrue(
            rel_error < 1E-2,
            msg="Actual error: {0}. Expected error: 1E-2".format(rel_error))
        self.assertTrue(
            len(res) < 10,
            msg="Needed {0} iterations to solve system.".format(len(res)) +
            "Expected not more than 9 iterations.")

    def test_laplace_potentials(self):
        """Test the Laplace potential operators."""

        grid = bempp.api.shapes.regular_sphere(4)

        space = bempp.api.function_space(grid, "P", 1)

        def dirichlet_fun(x):
            return np.exp(x[0]) * np.sin(x[1])

        def dirichlet_data(x, n, domain_index, res):
            res[0] = dirichlet_fun(x)

        def neumann_data(x, n, domain_index, res):
            grad = np.array([np.exp(x[0]) * np.sin(x[1]),
                             np.exp(x[0]) * np.cos(x[1]), 0])
            res[0] = np.dot(grad, n)

        dirichlet_grid_fun = bempp.api.GridFunction(space, fun=dirichlet_data)
        neumann_grid_fun = bempp.api.GridFunction(space, fun=neumann_data)

        #pylint: disable=no-member
        point = np.array([[0, 0.2, 0.3]]).T

        sl = bempp.api.operators.potential.laplace.single_layer(space, point)
        dl = bempp.api.operators.potential.laplace.double_layer(space, point)

        #pylint: disable=unsubscriptable-object
        actual = (sl * neumann_grid_fun - dl * dirichlet_grid_fun)[0, 0]
        expected = dirichlet_fun(point[:, 0])

        rel_error = np.abs(actual - expected) / np.abs(expected)
        self.assertTrue(
            rel_error < 1E-6,
            msg="Actual error: {0}. Expected error: 1E-6".format(rel_error))


if __name__ == "__main__":

    from unittest import main
    main()
