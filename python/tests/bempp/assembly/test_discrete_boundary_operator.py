import pytest
from bempp import grid_from_sphere
from bempp import function_space
from bempp.operators.boundary.laplace import single_layer as laplace_slp
from bempp.operators.boundary.helmholtz import single_layer as helmholtz_slp
from bempp.assembly.discrete_boundary_operator import ZeroDiscreteBoundaryOperator
from bempp.assembly.discrete_boundary_operator import BlockedDiscreteBoundaryOperator
from bempp import operators

import numpy as np

_eps = 1E-13


@pytest.fixture(scope='module')
def real_operator():
    grid = grid_from_sphere(3)
    space = function_space(grid,"DP",0)
    return laplace_slp(space,space,space).weak_form()

@pytest.fixture(scope='module')
def complex_operator():
    grid = grid_from_sphere(3)
    space = function_space(grid,"DP",0)
    return helmholtz_slp(space,space,space,1).weak_form()

@pytest.fixture(scope='module')
def zero_operator():
    return ZeroDiscreteBoundaryOperator(200,100)

class TestDiscreteBoundaryOperator(object):

    def test_correct_dtype(self,real_operator,complex_operator):

        assert real_operator.dtype=='float64'
        assert complex_operator.dtype=='complex128'

    def test_as_matrix_and_matvec_correct(self,real_operator,complex_operator):

        mat = np.eye(real_operator.shape[1],dtype='float64')
        expected_real = real_operator*mat
        assert expected_real.dtype=='float64'
        assert np.linalg.norm(real_operator.as_matrix()-expected_real,'fro')<_eps

        expected_complex = complex_operator*mat
        assert expected_complex.dtype=='complex128'
        assert np.linalg.norm(complex_operator.as_matrix()-expected_complex,'fro')<_eps

    def test_matvec_correct_dimension(self,real_operator):

        vec2d = np.eye(real_operator.shape[1],1,dtype='float64')
        vec1d = vec2d.ravel()

        res1d = real_operator*vec1d
        res2d = real_operator*vec2d

        assert res1d.ndim == 1
        assert res2d.ndim == 2

    def test_shape(self):

        grid = grid_from_sphere(2)
        space_const = function_space(grid,"DP",0)
        space_lin = function_space(grid,"P",1)
        op = laplace_slp(space_lin,space_lin,space_const).weak_form()
        assert op.shape==(space_const.global_dof_count,space_lin.global_dof_count)

class TestScaledDiscreteBoundaryOperator(object):

    @pytest.mark.parametrize('alpha',[2.0,2+1j])
    def test_correct_dtype(self,alpha,real_operator,complex_operator):

        real_scaled = alpha*real_operator
        complex_scaled = alpha*complex_operator

        res_type_1 = np.dtype(type(np.float64(1.0)*alpha))
        res_type_2 = np.dtype(type(np.complex128(1.0)*alpha))

        assert real_scaled.dtype==res_type_1
        assert complex_scaled.dtype==res_type_2

    @pytest.mark.parametrize('alpha',[2.0,2+1j])
    def test_as_matrix(self,alpha,real_operator,complex_operator):

        scaled_real_op = alpha*real_operator
        scaled_complex_op = alpha*complex_operator
        
        real_op_expected = alpha*(real_operator.as_matrix())
        real_op_actual = scaled_real_op.as_matrix()

        complex_op_expected = alpha*(complex_operator.as_matrix())
        complex_op_actual = scaled_complex_op.as_matrix()

        assert np.linalg.norm(real_op_expected-real_op_actual,'fro') < _eps
        assert np.linalg.norm(complex_op_expected-complex_op_actual,'fro') < _eps

    @pytest.mark.parametrize('alpha',[2.0,2+1j])
    def test_matvec(self,alpha,real_operator,complex_operator):

        scaled_real_op = alpha*real_operator
        scaled_complex_op = alpha*complex_operator

        vec = np.random.rand(scaled_real_op.shape[0])
        res_expected_real = alpha*(real_operator*vec)
        res_expected_complex = alpha*(complex_operator*vec)

        res_expected_actual_real = scaled_real_op*vec
        res_expected_actual_complex = scaled_complex_op*vec

        assert np.linalg.norm(res_expected_real-res_expected_actual_real)<_eps
        assert np.linalg.norm(res_expected_complex-res_expected_actual_complex)<_eps

    def test_shape(self,real_operator):

        alpha = 2.0
        scaled_real_op = alpha*real_operator

        assert scaled_real_op.shape==real_operator.shape


class TestSumDiscreteBoundaryOperator(object):

    def test_correct_dtype(self,real_operator,complex_operator):

        real_plus_real = real_operator+real_operator
        real_plus_complex = real_operator+complex_operator
        complex_plus_complex = complex_operator+complex_operator

        assert real_plus_real.dtype=='float64'
        assert real_plus_complex.dtype=='complex128'
        assert complex_plus_complex.dtype=='complex128'

    def test_as_matrix(self,real_operator,complex_operator):

        sum_operator = real_operator+complex_operator

        assert np.linalg.norm(sum_operator.as_matrix()-(
            real_operator.as_matrix()+complex_operator.as_matrix()))<_eps

    def test_matvec(self,real_operator,complex_operator):
        
        vec = np.random.rand(real_operator.shape[1])

        res_expected = real_operator*vec+complex_operator*vec
        res_actual = (real_operator+complex_operator)*vec

        assert np.linalg.norm(res_expected-res_actual)<_eps


    def test_shape(self,real_operator,complex_operator):

        sum_op = real_operator+complex_operator
        assert sum_op.shape==real_operator.shape

class TestZeroDiscreteBoundaryOperator(object):

    def test_correct_dtype(self,zero_operator):

        assert zero_operator.dtype=='float64'

    def test_as_matrix(self,zero_operator):

        mat = zero_operator.as_matrix()

        assert np.linalg.norm(mat)==0
        assert mat.shape==(200,100)

    def test_matvec(self,zero_operator):
        
        vec = np.random.rand(zero_operator.shape[1])

        res_expected = np.zeros((200,),dtype='float64')
        res_actual = zero_operator*vec

        assert res_actual.ndim==1
        assert np.linalg.norm(res_expected-res_actual)==0

    def test_shape(self,zero_operator):

        assert zero_operator.shape==(200,100)

class TestBlockedDiscreteBoundaryOperator(object):

    @pytest.fixture
    def row_dimensions(self):
        return np.array([20,512,40])

    @pytest.fixture
    def column_dimensions(self):
        return np.array([512,30])

    def test_initialization(self,row_dimensions,column_dimensions):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)

    def test_shape(self,row_dimensions,column_dimensions):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        assert op.shape == (572,542)

    def test_ndims(self,row_dimensions,column_dimensions):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        assert op.ndims == (3,2)

    def test_row_dimensions(self,row_dimensions,column_dimensions):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        assert np.linalg.norm(op.row_dimensions-row_dimensions)==0 

    def test_column_dimensions(self,row_dimensions,column_dimensions):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        assert np.linalg.norm(op.column_dimensions-column_dimensions)==0

    def test_set_operator(self,row_dimensions,column_dimensions,
            real_operator):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)

        op[1,0] = real_operator

    def test_dtype(self,row_dimensions,column_dimensions,
            real_operator,complex_operator):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)

        op[1,0] = real_operator

        assert op.dtype=='float64'

        op[1,0] = complex_operator

        assert op.dtype=='complex128'

    def test_matvec(self,row_dimensions,column_dimensions,
            real_operator):

        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)

        op[1,0] = real_operator

        vec = np.random.rand(op.shape[1],2)
        actual = op*vec

        expected = np.zeros((op.shape[0],vec.shape[1]),dtype=op.dtype)
        expected[20:20+512,:] = real_operator*vec[:512,:]
        assert np.linalg.norm(actual-expected)<_eps

    def test_as_matrix(self,row_dimensions,column_dimensions,real_operator):


        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        op[1,0] = real_operator

        expected = op*np.eye(op.shape[1],dtype=op.dtype)
        actual = op.as_matrix()

        assert np.linalg.norm(expected-actual)<_eps

    def test_scale(self,row_dimensions,column_dimensions,real_operator):

        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        op[1,0] = real_operator

        alpha = 1+2.0j
        res = alpha*op

        expected = alpha*(op.as_matrix())
        actual = res.as_matrix()

        assert np.linalg.norm(expected-actual)<_eps
        assert isinstance(res,BlockedDiscreteBoundaryOperator)

    def test_add_blocked_operator(self,row_dimensions,column_dimensions,
            real_operator):
        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        op[1,0] = real_operator

        res = op+op
        actual = res.as_matrix()
        expected = 2*op.as_matrix()

        assert isinstance(res,BlockedDiscreteBoundaryOperator)
        assert np.linalg.norm(actual-expected)<_eps

    def test_add_nonblocked_operator(self,row_dimensions,column_dimensions,
            real_operator):

        op = BlockedDiscreteBoundaryOperator(row_dimensions,column_dimensions)
        op[1,0] = real_operator

        res = op+ZeroDiscreteBoundaryOperator(op.shape[0],op.shape[1])
        actual = res.as_matrix()
        expected = op.as_matrix()

        assert np.linalg.norm(actual-expected)<_eps

class TestInverseSparseDiscreteBoundary(object):

    def test_inverse_of_square_matrix(self):
        
        from bempp.assembly.discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator
        grid = grid_from_sphere(3)
        space = function_space(grid,"P",1)
        sparse_op = operators.boundary.sparse.identity(space,space,space).weak_form()
        inverse_operator = InverseSparseDiscreteBoundaryOperator(sparse_op)
        actual = inverse_operator*sparse_op.sparse_operator.todense()
        expected = np.eye(*sparse_op.shape)
        assert np.linalg.norm(actual-expected)<_eps

    def test_inverse_of_thin_rectangular_matrix(self):

        from bempp.assembly.discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator
        grid = grid_from_sphere(3)
        trial_space= function_space(grid,"P",1)
        test_space = function_space(grid,"P",2)
        sparse_op = operators.boundary.sparse.identity(trial_space,trial_space,test_space).weak_form()
        print(sparse_op.shape)
        inverse_operator = InverseSparseDiscreteBoundaryOperator(sparse_op)
        actual = inverse_operator*sparse_op.sparse_operator.todense()
        expected = np.eye(*(inverse_operator.shape[0],sparse_op.shape[1]))
        assert np.linalg.norm(actual-expected)<_eps

    def test_inverse_of_thick_rectangular_matrix(self):

        from bempp.assembly.discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator
        grid = grid_from_sphere(3)
        trial_space= function_space(grid,"P",2)
        test_space = function_space(grid,"P",1)
        sparse_op = operators.boundary.sparse.identity(trial_space,trial_space,test_space).weak_form()
        print(sparse_op.shape)
        inverse_operator = InverseSparseDiscreteBoundaryOperator(sparse_op)
        actual = sparse_op*inverse_operator.as_matrix()
        expected = np.eye(*(sparse_op.shape[0],inverse_operator.shape[1]))
        assert np.linalg.norm(actual-expected)<_eps





        


