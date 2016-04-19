from unittest import TestCase
import unittest
import bempp.api
import numpy as np

class TestMaxwell(TestCase):

    def test_efie_unit_sphere_rt_functions(self):
	"""Solve an exterior EFIE problem un the unit sphere with Raviart-Thomas functions."""

	# This script solves the Maxwell equations in the region exterior to a bounded
	# object, with Dirichlet boundary conditions given by the exact solution
	# (satisfying the Silver-Mueller radiation conditions)
	#
	#     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
	#
	# where (r, theta, phi) are the radial, zenith angle and azimuthal spherical
	# coordinates in the system anchored at the point (0.1, 0.1, 0.1), h_1^{(1)}(r)
	# is the spherical Hankel function of the first kind and order 1 and \hat phi is
	# the unit vector oriented along d(\vec x)/d\phi.
	    
	k = 1
	source = 0.1

        grid = bempp.api.shapes.regular_sphere(4)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

       	def eval_dirichlet_data(point, normal, domain_index, result): 
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    scale = h1kr / r
	    field = [-y * scale, x * scale, 0.]
	    result[:] = np.cross(field, normal)

	def eval_exact_neumann_data(point, normal, domain_index, result):
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
			  np.exp(1j * kr) / (kr * kr * r))
	    xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
	    curl = [x * z * xy_factor,
		    y * z * xy_factor,
		    ((x*x + y*y + 2*z*z) * h1kr + r * (x*x + y*y) * h1kr_deriv) /
		    (r * r * r)]
	    result[:] = np.cross(curl, normal) / (1j * k)

	def eval_exact_solution(point):
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * exp(1j * kr) / (kr * kr)
	    scale = h1kr / r
	    return np.array([-y * scale, x * scale, 0.])		 

	space = bempp.api.function_space(grid, "RT", 0)

	efie = bempp.api.operators.boundary.maxwell.electric_field(space, space, space, k, parameters=parameters)
	mfie = bempp.api.operators.boundary.maxwell.magnetic_field(space, space, space, k, parameters=parameters)
	ident = bempp.api.operators.boundary.sparse.maxwell_identity(space, space, space, parameters=parameters)

	dirichlet_grid_fun = bempp.api.GridFunction(space, fun=eval_dirichlet_data)
	rhs_coeffs = -(.5 * ident.weak_form() + mfie.weak_form() ) * dirichlet_grid_fun.coefficients
	
        from scipy.linalg import solve
	sol_coefficients = solve(bempp.api.as_matrix(efie.weak_form()), rhs_coeffs)
        sol = bempp.api.GridFunction(space, coefficients=sol_coefficients)

	exact_solution = bempp.api.GridFunction(space, fun=eval_exact_neumann_data)
	rel_error = (sol-exact_solution).l2_norm() / exact_solution.l2_norm()
        self.assertTrue(rel_error  < 2E-2, msg="Actual error: {0}. Expected error: 2E-2".format(rel_error))

    def test_efie_unit_sphere_rwg_functions(self):
	"""Solve an exterior EFIE problem un the unit sphere with RWG functions."""

	# This script solves the Maxwell equations in the region exterior to a bounded
	# object, with Dirichlet boundary conditions given by the exact solution
	# (satisfying the Silver-Mueller radiation conditions)
	#
	#     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
	#
	# where (r, theta, phi) are the radial, zenith angle and azimuthal spherical
	# coordinates in the system anchored at the point (0.1, 0.1, 0.1), h_1^{(1)}(r)
	# is the spherical Hankel function of the first kind and order 1 and \hat phi is
	# the unit vector oriented along d(\vec x)/d\phi.
	    
	k = 1
	source = 0.1

        grid = bempp.api.shapes.regular_sphere(4)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'

       	def eval_dirichlet_data(point, normal, domain_index, result): 
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    scale = h1kr / r
	    field = [-y * scale, x * scale, 0.]
	    result[:] = np.cross(field, normal)

	def eval_exact_neumann_data(point, normal, domain_index, result):
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
			  np.exp(1j * kr) / (kr * kr * r))
	    xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
	    curl = [x * z * xy_factor,
		    y * z * xy_factor,
		    ((x*x + y*y + 2*z*z) * h1kr + r * (x*x + y*y) * h1kr_deriv) /
		    (r * r * r)]
	    result[:] = np.cross(curl, normal) / (1j * k)

	def eval_exact_solution(point):
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * exp(1j * kr) / (kr * kr)
	    scale = h1kr / r
	    return np.array([-y * scale, x * scale, 0.])		 

	space = bempp.api.function_space(grid, "RWG", 0)

	efie = bempp.api.operators.boundary.maxwell.electric_field(space, space, space, k, parameters=parameters)
	mfie = bempp.api.operators.boundary.maxwell.magnetic_field(space, space, space, k, parameters=parameters)
	ident = bempp.api.operators.boundary.sparse.maxwell_identity(space, space, space, parameters=parameters)

	dirichlet_grid_fun = bempp.api.GridFunction(space, fun=eval_dirichlet_data)
	rhs_coeffs = -(.5 * ident.weak_form() + mfie.weak_form() ) * dirichlet_grid_fun.coefficients
	
        from scipy.linalg import solve
	sol_coefficients = solve(bempp.api.as_matrix(efie.weak_form()), rhs_coeffs)
        sol = bempp.api.GridFunction(space, coefficients=sol_coefficients)

	exact_solution = bempp.api.GridFunction(space, fun=eval_exact_neumann_data)
	rel_error = (sol-exact_solution).l2_norm() / exact_solution.l2_norm()
        self.assertTrue(rel_error  < 2E-2, msg="Actual error: {0}. Expected error: 2E-2".format(rel_error))


    def test_maxwell_potential_operators(self):
	"""Solve an exterior EFIE problem un the unit sphere with RWG functions."""

	# This script solves the Maxwell equations in the region exterior to a bounded
	# object, with Dirichlet boundary conditions given by the exact solution
	# (satisfying the Silver-Mueller radiation conditions)
	#
	#     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
	#
	# where (r, theta, phi) are the radial, zenith angle and azimuthal spherical
	# coordinates in the system anchored at the point (0.1, 0.1, 0.1), h_1^{(1)}(r)
	# is the spherical Hankel function of the first kind and order 1 and \hat phi is
	# the unit vector oriented along d(\vec x)/d\phi.
	    
	k = 1
	source = 0.1

        grid = bempp.api.shapes.regular_sphere(5)
        parameters = bempp.api.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'dense'
        parameters.assembly.potential_operator_assembly_type = 'dense'

       	def eval_dirichlet_data(point, normal, domain_index, result): 
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    scale = h1kr / r
	    field = [-y * scale, x * scale, 0.]
	    result[:] = np.cross(field, normal)

	def eval_exact_neumann_data(point, normal, domain_index, result):
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
			  np.exp(1j * kr) / (kr * kr * r))
	    xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
	    curl = [x * z * xy_factor,
		    y * z * xy_factor,
		    ((x*x + y*y + 2*z*z) * h1kr + r * (x*x + y*y) * h1kr_deriv) /
		    (r * r * r)]
	    result[:] = np.cross(curl, normal) / (1j * k)

	def eval_exact_solution(point):
	    x, y, z = point - source
	    r = np.sqrt(x**2 + y**2 + z**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    scale = h1kr / r
	    return np.array([-y * scale, x * scale, 0.])		 


        space = bempp.api.function_space(grid, "RWG", 0)

	dirichlet_data = bempp.api.GridFunction(space, fun=eval_dirichlet_data)
        neumann_data = bempp.api.GridFunction(space, fun=eval_exact_neumann_data)

        eval_point = np.array([[3, 2, 1]]).T

        efie_pot = bempp.api.operators.potential.maxwell.electric_field(space, eval_point, k,
                parameters=parameters)
        mfie_pot = bempp.api.operators.potential.maxwell.magnetic_field(space, eval_point, k,
                parameters=parameters)
        
        expected = eval_exact_solution(eval_point[:, 0])
        actual = (-efie_pot * neumann_data - mfie_pot * dirichlet_data)[:, 0]
        rel_error = np.linalg.norm(expected-actual) / np.linalg.norm(actual)

        self.assertTrue(rel_error  < 1E-3, msg="Actual error: {0}. Expected error: 1E-3".format(rel_error))
    
    def test_electric_far_field(self):
        """Test the electric far field operator."""

	# This example computes the far-field pattern for the solution of
	# $$
	# \text{curl}~\text{curl} E - k^2 E = 0
	# $$
	# with boundary data given from an analytic Maxwell solution in the exterior of the sphere as
	# $$
	# E\times n = h_1^{(1)}(kr)e_{\phi}\times n
	# $$
	# Here, $e_{\phi}$ is the unit vector along the $\phi$ derivative $\frac{dx}{d\phi}$ in spherical $(r,\theta,\phi)$ coordinates.
	#
	# The far-field pattern is given analytically by $-e_{\phi}$.

	bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'dense'
        bempp.api.global_parameters.assembly.potential_operator_assembly_type = 'dense'
	grid = bempp.api.shapes.regular_sphere(4)

	space = bempp.api.function_space(grid,'RWG',0)

	k = 1

	def dirichlet_data(x, n, domain_index, res):
	    r = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
	    kr = k * r
	    h1kr = (-1j - kr) * np.exp(1j * kr) / (kr * kr)
	    field = h1kr * np.array([-x[1] / r, x[0] / r, 0.])
	    res[:] = np.cross(field, n)

	grid_fun = bempp.api.GridFunction(space, fun=dirichlet_data)

	ident = bempp.api.as_matrix(
                bempp.api.operators.boundary.sparse.maxwell_identity(space, space, space).weak_form())
	efie = bempp.api.as_matrix(
                bempp.api.operators.boundary.maxwell.electric_field(space, space, space, k).weak_form())

	from scipy.linalg import solve
	sol_coeffs = solve(efie, ident * grid_fun.coefficients)
	sol = bempp.api.GridFunction(space, coefficients=sol_coeffs)

	from bempp.api.operators.far_field.maxwell import electric_field as electric_far_field
	npoints = 100
	theta = np.linspace(0, 2 * np.pi, npoints)
	points = np.vstack([np.cos(theta), np.sin(theta), np.zeros(100,dtype='float64')])

	far_field = electric_far_field(space, points, k) * sol
	exact = np.vstack([points[1], -points[0], np.zeros(100)])

	rel_error = np.linalg.norm(far_field-exact)/np.linalg.norm(exact)
        self.assertTrue(rel_error  < 1E-5, msg="Actual error: {0}. Expected error: 1E-5".format(rel_error))

if __name__ == "__main__":

    from unittest import main
    main()


