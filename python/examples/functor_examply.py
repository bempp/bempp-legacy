import sys
import numpy
sys.path.append("..")
from bempp import surfaceNormalIndependentFunctor,testFunctor,SurfaceNormalIndependentFunctor_float64



def fun(point):
    return point[0]

if __name__=="__main__":
    f=surfaceNormalIndependentFunctor(fun)
    res=numpy.array([0],dtype='float64')
    f.evaluate([1.,2.,3.],res)
    print testFunctor(f)
    print res
