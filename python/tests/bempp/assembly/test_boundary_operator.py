import pytest

from bempp import grid_from_sphere
from bempp import function_space
from bempp.operators.boundary.laplace import single_layer as laplace_slp
from bempp.operators.boundary.helmholtz import single_layer as helmholtz_slp
from bempp.operators.boundary.sparse import identity
from bempp.assembly.boundary_operator import _SumBoundaryOperator
from bempp.assembly.boundary_operator import _ScaledBoundaryOperator
from bempp.assembly.boundary_operator import BlockedBoundaryOperator
from bempp import GridFunction

import bempp


import numpy as np

_eps = 1E-13

@pytest.fixture(scope='module')
def space():
    grid = grid_from_sphere(3)
    space = function_space(grid,"DP",0)
    return space

@pytest.fixture(scope='module')
def real_operator(space):
    return laplace_slp(space,space,space)

@pytest.fixture(scope='module')
def complex_operator(space):
    return helmholtz_slp(space,space,space,1)

@pytest.fixture(scope='module')
def zero_operator(space):
    return ZeroBoundaryOperator(space,space,space)

@pytest.fixture(scope='module')
def sparse_operator(space):
    return identity(space,space,space)

class TestBoundaryOperator(object):

    def test_weak_form(self,real_operator,complex_operator,space):

        real_weak = real_operator.weak_form()
        complex_weak = complex_operator.weak_form()

        assert real_weak.shape==(space.global_dof_count,space.global_dof_count)
        assert complex_weak.shape==(space.global_dof_count,space.global_dof_count)


class TestSumBoundaryOperator(object):

    def test_add_boundary_operators(self,real_operator,complex_operator):

        res_real = real_operator+real_operator
        res_complex = real_operator+complex_operator

        assert res_real.result_type == 'float64'
        assert res_complex.result_type == 'complex128'

        expected_real = real_operator.weak_form()+real_operator.weak_form()
        expected_complex = real_operator.weak_form()+complex_operator.weak_form()

        actual_real = res_real.weak_form()
        actual_complex = res_complex.weak_form()

        assert isinstance(res_real,_SumBoundaryOperator)
        assert np.linalg.norm(expected_real.as_matrix()-actual_real.as_matrix())<_eps
        assert np.linalg.norm(expected_complex.as_matrix()-actual_complex.as_matrix())<_eps


class TestScaledBoundaryOperator:

    def test_scale_boundary_operator(self,real_operator):

        alpha_real = 2.0
        alpha_complex = 1.0+2.0j

        res_real = alpha_real*real_operator
        res_complex = alpha_complex*real_operator

        assert res_real.result_type == 'float64'
        assert res_complex.result_type == 'complex128'

        assert isinstance(res_real,_ScaledBoundaryOperator)

        actual_real = res_real.weak_form().as_matrix()
        expected_real = alpha_real*real_operator.weak_form().as_matrix()

        actual_complex = res_complex.weak_form().as_matrix()
        expected_complex = alpha_complex*real_operator.weak_form().as_matrix()

        assert np.linalg.norm(actual_real-expected_real)<_eps
        assert np.linalg.norm(actual_complex-expected_complex)<_eps

class TestBlockedBoundaryOperator:

    def test_initialization(self):

        blocked_operator = BlockedBoundaryOperator(3,2)

    def test_assign_operators(self,real_operator,complex_operator):

        blocked_operator = BlockedBoundaryOperator(3,2)

        blocked_operator[0,0] = real_operator
        assert blocked_operator.domain_spaces[0] is not None
        assert blocked_operator.range_spaces[0] is not None
        assert blocked_operator.dual_to_range_spaces[0] is not None
        assert blocked_operator.domain_spaces[1] is None

    def test_result_type(self,real_operator,complex_operator):

        blocked_operator = BlockedBoundaryOperator(3,2)
        blocked_operator[0,0] = real_operator
        assert blocked_operator.result_type=='float64'

        blocked_operator[1,1] = complex_operator
        assert blocked_operator.result_type=='complex128'

    def test_basis_type(self):

        blocked_operator = BlockedBoundaryOperator(3,2)
        assert blocked_operator.basis_type=='float64'

    def test_weak_form(self,real_operator,complex_operator):

        blocked_operator = BlockedBoundaryOperator(3,2)
        blocked_operator[0,0] = real_operator
        with pytest.raises(ValueError):
            weak_form = blocked_operator.weak_form()
        blocked_operator[1,1] = complex_operator
        blocked_operator[2,0] = real_operator
        weak_form = blocked_operator.weak_form()

        assert weak_form.ndims == (3,2)

    def test_apply_grid_function(self,real_operator,complex_operator,space):

        blocked_operator = BlockedBoundaryOperator(2,2)
        blocked_operator[0,0] = real_operator
        blocked_operator[1,0] = complex_operator
        blocked_operator[1,1] = real_operator

        g = GridFunction(space=space,
                coefficients = np.ones(space.global_dof_count))
        g_vec = [g,2*g]
        res = blocked_operator*g_vec
        actual = res[1].coefficients
        expected = (2*real_operator*g+complex_operator*g).coefficients

        assert np.linalg.norm(actual-expected)<_eps

    def test_add_blocked_boundary_operator(self,real_operator,complex_operator):

        blocked_operator = BlockedBoundaryOperator(2,2)
        blocked_operator[0,0] = real_operator
        blocked_operator[1,0] = complex_operator
        blocked_operator[1,1] = real_operator

        result = blocked_operator+blocked_operator

        expected = 2*blocked_operator.weak_form().as_matrix()
        actual = result.weak_form().as_matrix()

        assert np.linalg.norm(expected-actual)<_eps

class TestProductBoundaryOperator(object):

    def test_product_weak_form_agrees_with_product_of_weak_forms(self):

        grid = bempp.grid_from_sphere(3)
        space = bempp.function_space(grid,"P",1)
        slp = bempp.operators.boundary.laplace.single_layer(space,space,space)
        hyp = bempp.operators.boundary.laplace.single_layer(space,space,space)
        op = slp*hyp
        actual = op.weak_form().as_matrix()
        expected = hyp.weak_form().as_matrix()*slp.strong_form().as_matrix()
        assert np.linalg.norm(actual-expected)<_eps

    def test_apply_grid_function(self):

        grid = bempp.grid_from_sphere(3)
        space = bempp.function_space(grid,"P",1)
        slp = bempp.operators.boundary.laplace.single_layer(space,space,space)
        hyp = bempp.operators.boundary.laplace.single_layer(space,space,space)
        op = slp*hyp

        g = GridFunction(slp.domain,coefficients=np.ones(slp.domain.global_dof_count))
        expected = hyp*(slp*g)
        actual = op*g

        assert np.linalg.norm(expected.coefficients-actual.coefficients)<_eps



