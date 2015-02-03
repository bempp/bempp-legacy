import pytest

from bempp import grid_from_sphere
from bempp import function_space
from bempp.operators.boundary.laplace import single_layer as laplace_slp
from bempp.operators.boundary.helmholtz import single_layer as helmholtz_slp
from bempp.operators.boundary.sparse import identity
from bempp.assembly.boundary_operator import _SumBoundaryOperator
from bempp.assembly.boundary_operator import _ScaledBoundaryOperator


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

    def test_label(self,real_operator):

        res = real_operator+real_operator
        label = res.label
        assert isinstance(label,str)

class TestScaledBoundaryOperator(object):

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

    def test_label(self,real_operator):

        alpha = 2.0

        res = alpha*real_operator

        label = res.label

        assert isinstance(label,str)


