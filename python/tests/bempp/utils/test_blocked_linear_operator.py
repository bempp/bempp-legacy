from py.test import mark, fixture
from bempp.utils import BlockedDiscreteLinearOperator, BlockedLinearOperator


class TestBlockedLinearOperator(object):

    @fixture
    def operator(self):

        from bempp import grid_from_sphere, function_space
        from bempp.operators.boundary.sparse import identity
        grid = grid_from_sphere(2)
        space = function_space(grid,"DP",0)
        return identity(space,space,space)

    def test_instantiate(self,operator):

        blocked_operator = BlockedLinearOperator(2,2)
        blocked_operator[0,0] = operator
        blocked_operator[1,1] = operator

        assert blocked_operator._rows[0]==True
        assert blocked_operator._rows[1]==True
        assert blocked_operator._cols[0]==True
        assert blocked_operator._cols[1]==True

    def test_weak_form(self,operator):

        blocked_operator = BlockedLinearOperator(2,2)
        blocked_operator[0,0] = operator
        blocked_operator[1,1] = operator

        weak_form = blocked_operator.weak_form()

        assert weak_form.ndims == (2,2)
        assert weak_form[0,0] is not None
        assert weak_form[1,1] is not None
        assert weak_form[0,1] is None
        assert weak_form[1,0] is None


class TestBlockedDiscreteLinearOpeartor(object):

    def real_operator(self,m,n):
        from scipy import sparse
        import numpy as np
        return sparse.rand(m,n,0.1,dtype=np.float64)

    def complex_operator(self,m,n):
        from scipy import sparse
        import numpy as np
        return (1+1j)*sparse.rand(m,n,0.1,dtype=np.float64)

    def instantiate_real_operator(self):

        op1 = self.real_operator(100,100)
        op2 = self.real_operator(100,200)
        op3 = self.real_operator(100,200)

        blocked_operator = BlockedDiscreteLinearOperator(2,2)
        blocked_operator[0,0] = op1
        blocked_operator[0,1] = op2
        blocked_operator[1,1] = op3

        
        return blocked_operator

    def instantiate_complex_operator(self):

        op1 = self.complex_operator(100,100)
        op2 = self.real_operator(100,200)
        op3 = self.real_operator(100,200)

        blocked_operator = BlockedDiscreteLinearOperator(2,2)
        blocked_operator[0,0] = op1
        blocked_operator[0,1] = op2
        blocked_operator[1,1] = op3
        
        return blocked_operator

    def test_shape(self):

        blocked_operator = self.instantiate_real_operator()

        assert blocked_operator.shape == (200,300)

    def test_dtype(self):

        real_blocked_operator = self.instantiate_real_operator()
        complex_blocked_operator = self.instantiate_complex_operator()
        
        assert real_blocked_operator.dtype == 'float64'
        assert complex_blocked_operator.dtype == 'complex128'

    def test_matvec(self):

        real_blocked_operator = self.instantiate_real_operator()
        complex_blocked_operator = self.instantiate_complex_operator()

        import numpy as np

        cols = real_blocked_operator.shape[1] 
        rows = real_blocked_operator.shape[0]

        x_real = np.random.rand(cols)
        x_complex = np.random.rand(cols)+1j*np.random.rand(cols)

        expected_real = np.zeros(rows,dtype=np.float64)
        expected_real[:100] = real_blocked_operator[0,0]*x_real[:100]+real_blocked_operator[0,1]*x_real[100:]
        expected_real[100:] = real_blocked_operator[1,1]*x_real[100:]

        actual_real = real_blocked_operator*x_real

        assert np.linalg.norm(actual_real-expected_real)<1E-14

        expected_complex = np.zeros(rows,dtype=np.complex128)
        expected_complex[:100] = real_blocked_operator[0,0]*x_complex[:100]+real_blocked_operator[0,1]*x_complex[100:]
        expected_complex[100:] = real_blocked_operator[1,1]*x_complex[100:]

        actual_complex = real_blocked_operator*x_complex
        assert np.linalg.norm(actual_complex-expected_complex,np.inf)<1E-14

    def test_scalar_mult(self):


        real_blocked_operator = self.instantiate_real_operator()
        alpha = 1+1j

        import numpy as np

        cols = real_blocked_operator.shape[1] 

        x = np.random.rand(cols)

        expected = alpha*(real_blocked_operator*x)
        actual = (alpha*real_blocked_operator)*x

        assert np.linalg.norm(actual-expected,np.inf)<1E-14







        

        

