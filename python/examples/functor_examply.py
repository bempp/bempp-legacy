import sys
import numpy
sys.path.append("..")
from bempp import PythonSurfaceNormalIndependentFunctor_float64,testFunctor



def fun(point):
    return point[0]

if __name__=="__main__":
    f=PythonSurfaceNormalIndependentFunctor_float64(fun,3,1)
    res=numpy.array([0],dtype='float64')
    f.evaluate([1.,2.,3.],res)
    print testFunctor(f)
