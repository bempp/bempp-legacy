import pytest
from bempp import grid_from_sphere
from bempp import function_space
from bempp import GridFunction
import numpy as np

_eps = 1E-13

class TestGridFunction(object):

    @pytest.fixture
    def grid(self):
        return grid_from_sphere(3)

    @pytest.fixture
    def space(self,grid):
        return function_space(grid,"DP",0)

    @pytest.fixture
    def dual_space(self,grid):
        return function_space(grid,"P",1)

    def test_instantiate_from_coefficients(self,space):
        coefficients = np.random.rand(space.global_dof_count)
        fun = GridFunction(space,coefficients=coefficients)
        assert np.linalg.norm(coefficients-fun.coefficients)==0

    def test_instantiate_from_projections(self,space,dual_space):
        projections = np.random.rand(dual_space.global_dof_count)
        fun = GridFunction(space,dual_space=dual_space,projections=projections)
        assert np.linalg.norm(projections-fun.projections(dual_space))<_eps

    def test_instantiate_from_python_function(self,space,dual_space):

        def py_fun(x,n,domain_index,res):
            res[0] = 1
        coefficients = np.ones(space.global_dof_count)
        fun = GridFunction(space,dual_space=space,fun=py_fun)
        assert np.linalg.norm(coefficients-fun.coefficients)<_eps

    def test_instantiate_from_complex_python_function(self,space,dual_space):

        def py_fun(x,n,domain_index,res):
            res[0] = 1j
        coefficients = 1j*np.ones(space.global_dof_count)
        fun = GridFunction(space,dual_space=space,fun=py_fun,complex_data=True)
        assert np.linalg.norm(coefficients-fun.coefficients)<_eps

    def test_add_grid_functions(self,space):

        coefficients = np.random.rand(space.global_dof_count)
        fun = GridFunction(space,coefficients=coefficients)
        res = fun+fun
        expected = 2*fun.coefficients
        actual = res.coefficients

        assert np.linalg.norm(expected-actual)<_eps






