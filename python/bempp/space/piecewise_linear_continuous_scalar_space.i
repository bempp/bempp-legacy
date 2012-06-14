%{
  #include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"


%include "space/piecewise_linear_continuous_scalar_space.hpp"

namespace Bempp
{

  %template(PiecewiseLinearContinuousScalarSpace_float64)     PiecewiseLinearContinuousScalarSpace<double>;
  %template(PiecewiseLinearContinuousScalarSpace_float32)     PiecewiseLinearContinuousScalarSpace<float>;
  %template(PiecewiseLinearContinuousScalarSpace_complex64)   PiecewiseLinearContinuousScalarSpace<std::complex<float> >;
  %template(PiecewiseLinearContinuousScalarSpace_complex128)  PiecewiseLinearContinuousScalarSpace<std::complex<double> >;

}

%pythoncode %{

class PiecewiseLinearContinuousScalarSpace(Template1,ScalarSpace):
    """Space of piecewise linear scalar functions"""
    def __init__(self,dtype1,*args,**kwargs):
        super(PiecewiseLinearContinuousScalarSpace,self).__init__('PiecewiseLinearContinuousScalarSpace',dtype1,*args,**kwargs)

	  %}

