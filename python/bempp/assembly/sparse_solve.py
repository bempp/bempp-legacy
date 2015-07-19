import numpy as np

class SparseSolve(object):

    def __init__(self,op):

        from scipy.sparse import csc_matrix
        from scipy.sparse.linalg import splu

        self._op = op
        if not isinstance(op,csc_matrix):
            raise ValueError("op must be of type scipy.sparse.csc.csc_matrix")

        if op.shape[0]==op.shape[1]:
            # Square matrix case
            solver = splu(op)
            self._solve_fun = lambda x: solver.solve(x)
        elif op.shape[0]>op.shape[1]:
            # Thin matrix case

            op_hermitian = op.conjugate().transpose()
            solver = splu((op_hermitian*op).tocsc())
            self._solve_fun = lambda x: solver.solve(op_hermitian*x)
        else:
            # Thick matrix case
            
            op_hermitian = op.conjugate().transpose()
            solver = splu((op*op_hermitian).tocsc())
            self._solve_fun = lambda x: op_hermitian*solver.solve(x)

    def solve(self,x):

        from bempp.utils.data_types import combined_type

        # Need self._op.shape[1] as that's right dimension for the inverse
        res = np.zeros((self._op.shape[1],x.shape[1]),combined_type(x.dtype,self._op.dtype))
        for i in range(x.shape[1]):
            # Ugly hack as Cython passes x as compiled ndarray
            # Need standard ndarray as otherwise the squeeze
            # function does not return 1-d array.
            vec = np.array(x[:,i]).squeeze()
            res[:,i] = self._solve_fun(vec)

        return res

        
