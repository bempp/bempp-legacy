import pytest
import bempp
import numpy as np

# Make sure that all quadrature orders are equal to compare dense and hmat results.

class TestHMatBoundaryOperator(object):

    @pytest.fixture
    def hmat_dense_compression_parameters(self):

        parameters = bempp.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'hmat'
        parameters.hmat.compression_algorithm = 'dense'

        parameters.quadrature.near.double_order = 4
        parameters.quadrature.medium.double_order = 4
        parameters.quadrature.far.double_order = 4

        return parameters

    @pytest.fixture
    def dense_operator_parameters(self):

        parameters = bempp.common.global_parameters()
        parameters.quadrature.near.double_order = 4
        parameters.quadrature.medium.double_order = 4
        parameters.quadrature.far.double_order = 4

        return parameters

    @pytest.fixture
    def grid(self):

        return bempp.grid_from_sphere(3)

    @pytest.fixture
    def const_space(self,grid):

        return bempp.function_space(grid,"DP",0)

    @pytest.fixture
    def random_vector(self,const_space):
        
        return np.random.randn(const_space.global_dof_count)


    @pytest.fixture
    def hmat_laplace_single_layer_operator(self,const_space, 
            hmat_dense_compression_parameters):

        from bempp import operators

        return operators.boundary.laplace.single_layer(
                const_space,const_space,const_space,
                parameters=hmat_dense_compression_parameters)

    @pytest.fixture
    def dense_laplace_single_layer_operator(self,const_space, 
            dense_operator_parameters):

        from bempp import operators

        return operators.boundary.laplace.single_layer(
                const_space,const_space,const_space,
                parameters=dense_operator_parameters)

    @pytest.fixture
    def hmat_laplace_single_layer_weak_form(self,
            hmat_laplace_single_layer_operator):

        return hmat_laplace_single_layer_operator.weak_form()

    @pytest.fixture
    def dense_laplace_single_layer_weak_form(self,
            dense_laplace_single_layer_operator):

        return dense_laplace_single_layer_operator.weak_form()


    def test_laplace_single_layer_as_matrix(self,
            hmat_laplace_single_layer_weak_form,
            dense_laplace_single_layer_weak_form):

        from bempp import operators

        expected = dense_laplace_single_layer_weak_form.as_matrix()
        actual = hmat_laplace_single_layer_weak_form.as_matrix()

        assert np.linalg.norm(expected-actual,np.Inf)<1E-15

    def test_laplace_single_layer_vec_mult(self,
            hmat_laplace_single_layer_weak_form,
            dense_laplace_single_layer_weak_form,
            random_vector):

        expected = dense_laplace_single_layer_weak_form*random_vector
        actual = hmat_laplace_single_layer_weak_form*random_vector

        assert np.linalg.norm(expected-actual,np.Inf)<1E-15











        





        
        


