%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

//namespace Bempp
//{
//  template <typename BasisFunctionType> class PiecewiseConstantScalarSpace : 
//}


%include "space/piecewise_constant_scalar_space.hpp"

namespace Bempp
{

  %template(PiecewiseConstantScalarSpace_float64)     PiecewiseConstantScalarSpace<double>;
  %template(PiecewiseConstantScalarSpace_float32)     PiecewiseConstantScalarSpace<float>;
  %template(PiecewiseConstantScalarSpace_complex64)   PiecewiseConstantScalarSpace<std::complex<float> >;
  %template(PiecewiseConstantScalarSpace_complex128)  PiecewiseConstantScalarSpace<std::complex<double> >;

}

%pythoncode %{

class PiecewiseConstantScalarSpace(Template1,ScalarSpace):
    """Space of piecewise constant scalar functions"""
    def __init__(self,dtype1,*args,**kwargs):
        super(PiecewiseConstantScalarSpace,self).__init__('PiecewiseConstantScalarSpace',dtype1,*args,**kwargs)

	  %}

