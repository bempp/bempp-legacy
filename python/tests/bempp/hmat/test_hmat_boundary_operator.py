import pytest
import bempp
import numpy as np

# Make sure that all quadrature orders are equal to compare dense and hmat results.

ACA_ACCURACY = 1E-3

class TestHMatBoundaryOperator(object):


    @pytest.fixture(params=['dense','aca'])
    def hmat_compression_parameters(self,request):

        parameters = bempp.common.global_parameters()
        parameters.assembly.boundary_operator_assembly_type = 'hmat'
        parameters.hmat.compression_algorithm = request.param
        parameters.hmat.eps = ACA_ACCURACY

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

        return bempp.grid_from_sphere(4)

    @pytest.fixture
    def const_space(self,grid):

        return bempp.function_space(grid,"DP",0)

    @pytest.fixture
    def random_vector(self,const_space):
        
        return np.random.randn(const_space.global_dof_count)


    @pytest.fixture
    def hmat_laplace_single_layer_operator(self,const_space, 
            hmat_compression_parameters):

        from bempp import operators

        return operators.boundary.laplace.single_layer(
                const_space,const_space,const_space,
                parameters=hmat_compression_parameters)

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
            dense_laplace_single_layer_weak_form,
            hmat_compression_parameters):
        pass

        #from bempp import operators

        #expected = dense_laplace_single_layer_weak_form.as_matrix()
        #actual = hmat_laplace_single_layer_weak_form.as_matrix()

        #if hmat_compression_parameters.hmat.compression_algorithm == 'aca':
        #    assert np.linalg.norm(expected-actual)/np.linalg.norm(expected)<ACA_ACCURACY
        #elif hmat_compression_parameters.hmat.compression_algorithm == 'dense':
        #    assert np.linalg.norm(expected-actual)/np.linalg.norm(expected)<1E-15

    def test_laplace_single_layer_vec_mult(self,
            hmat_laplace_single_layer_weak_form,
            dense_laplace_single_layer_weak_form,
            random_vector,
            hmat_compression_parameters):


        expected = dense_laplace_single_layer_weak_form*random_vector
        actual = hmat_laplace_single_layer_weak_form*random_vector

        if hmat_compression_parameters.hmat.compression_algorithm == 'aca':
            assert np.linalg.norm(expected-actual)/np.linalg.norm(expected)< ACA_ACCURACY
        elif hmat_compression_parameters.hmat.compression_algorithm == 'dense':
            assert np.linalg.norm(expected-actual)/np.linalg.norm(expected)<1E-15










        





        
        


