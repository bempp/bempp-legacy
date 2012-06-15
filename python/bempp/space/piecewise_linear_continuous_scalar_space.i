%{
  #include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"



%feature("pythonappend") Bempp::PiecewiseLinearContinuousScalarSpace<float>::PiecewiseLinearContinuousScalarSpace(Grid& grid)%{
self._BasisFunctionType='float32'
  %}


%feature("pythonappend") Bempp::PiecewiseLinearContinuousScalarSpace<double>::PiecewiseLinearContinuousScalarSpace(Grid& grid)%{
self._BasisFunctionType='float64'
  %}

%feature("pythonappend") Bempp::PiecewiseLinearContinuousScalarSpace<std::complex<float> >::PiecewiseLinearContinuousScalarSpace(Grid& grid)%{
self._BasisFunctionType='complex64'
  %}

%feature("pythonappend") Bempp::PiecewiseLinearContinuousScalarSpace<std::complex<double> >::PiecewiseLinearContinuousScalarSpace(Grid& grid)%{
self._BasisFunctionType='complex128'
  %}



%include "space/piecewise_linear_continuous_scalar_space.hpp"





namespace Bempp
{

  %template(PiecewiseLinearContinuousScalarSpace_float64)     PiecewiseLinearContinuousScalarSpace<double>;
  %template(PiecewiseLinearContinuousScalarSpace_float32)     PiecewiseLinearContinuousScalarSpace<float>;
  %template(PiecewiseLinearContinuousScalarSpace_complex64)   PiecewiseLinearContinuousScalarSpace<std::complex<float> >;
  %template(PiecewiseLinearContinuousScalarSpace_complex128)  PiecewiseLinearContinuousScalarSpace<std::complex<double> >;

}

%pythoncode %{


def PiecewiseLinearContinuousScalarSpace(grid,BasisFunctionType='float64'):
   """Create a space if piecewise linear constant functions."""
   name='PiecewiseLinearContinuousScalarSpace'+'_'+checkType(BasisFunctionType)
   return globals()[name](grid)

	  %}

