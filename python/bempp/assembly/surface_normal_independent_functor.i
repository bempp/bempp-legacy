%{
#include "assembly/surface_normal_independent_functor.hpp"
  %}



%define SURFACE_NORMAL_INDEPENDENT_FUNCTOR_REAL(TYPE,NPYNAME)
%apply arma::Col<TYPE>& INPLACE_COL { arma::Col<TYPE>& result_ };
%apply const arma::Col<TYPE>& IN_COL { const arma::Col<TYPE>& point };


%feature("director") SurfaceNormalIndependentFunctor<TYPE>;
%include "assembly/surface_normal_independent_functor.hpp"


namespace Bempp{
  %template(SurfaceNormalIndependentFunctor_## NPYNAME) SurfaceNormalIndependentFunctor<TYPE>;
    }

%clear arma::Col<TYPE>& result_;
%clear const arma::Col<TYPE>& point;

%enddef // SURFACE_NORMAL_INDEPENDENT_FUNCTOR_REAL

%define SURFACE_NORMAL_INDEPENDENT_FUNCTOR_COMPLEX(TYPE,NPYNAME)
%apply arma::Col<std::complex<TYPE> >& ARGOUT_COL { arma::Col<std::complex<TYPE> >& result_ };
%apply const arma::Col<TYPE>& IN_COL { const arma::Col<TYPE>& point };


%feature("director") SurfaceNormalIndependentFunctor<std::complex<TYPE> >;
%include "assembly/surface_normal_independent_functor.hpp"


namespace Bempp{
  %template(SurfaceNormalIndependentFunctor_## NPYNAME) SurfaceNormalIndependentFunctor<std::complex<TYPE> >;
    }

%clear arma::Col<std::complex<TYPE> >& result_;
%clear const arma::Col<std::complex<TYPE> >& point;

%enddef // SURFACE_NORMAL_INDEPENDENT_FUNCTOR_REAL


SURFACE_NORMAL_INDEPENDENT_FUNCTOR_REAL(float,float32)
SURFACE_NORMAL_INDEPENDENT_FUNCTOR_REAL(double,float64)
SURFACE_NORMAL_INDEPENDENT_FUNCTOR_COMPLEX(float,complex32)
SURFACE_NORMAL_INDEPENDENT_FUNCTOR_COMPLEX(double,complex64)

  %pythoncode %{

def surfaceNormalIndependentFunctor(fun,valueType='float64',argumentDimension=3,resultDimension=1):

    dtype=checkType(valueType)
    base=None
    if dtype=='float32':
        base=SurfaceNormalIndependentFunctor_float32
    elif dtype=='float64':
        base=SurfaceNormalIndependentFunctor_float64
    elif dtype=='complex64':
        base=SurfaceNormalIndependentFunctor_complex64
    elif dtype=='complex128':
        base=SurfaceNormalIndependentFunctor_complex128
    class SurfaceNormalIndependentFunctor(base):
        def __init__(self,fun,argumentDimension,resultDimension):
            self.fun=fun
            self.argDim=argumentDimension
            self.resultDim=resultDimension
	    self._valueType=dtype
	    super(SurfaceNormalIndependentFunctor,self).__init__()
        def argumentDimension(self):
            return self.argDim
        def resultDimension(self):
            return self.resultDim
        def evaluate(self,point,result):
            result[:]=self.fun(point) 
    return SurfaceNormalIndependentFunctor(fun,argumentDimension,resultDimension)

	      %}









