%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

//namespace Bempp
//{
//  template <typename BasisFunctionType> class PiecewiseConstantScalarSpace : 
//}


%feature("pythonappend") Bempp::PiecewiseConstantScalarSpace<float>::PiecewiseConstantScalarSpace(Grid& grid)%{
self._BasisFunctionType='float32'
  %}


%feature("pythonappend") Bempp::PiecewiseConstantScalarSpace<double>::PiecewiseConstantScalarSpace(Grid& grid)%{
self._BasisFunctionType='float64'
  %}

%feature("pythonappend") Bempp::PiecewiseConstantScalarSpace<std::complex<float> >::PiecewiseConstantScalarSpace(Grid& grid)%{
self._BasisFunctionType='complex64'
  %}

%feature("pythonappend") Bempp::PiecewiseConstantScalarSpace<std::complex<double> >::PiecewiseConstantScalarSpace(Grid& grid)%{
self._BasisFunctionType='complex128'
  %}

%include "space/piecewise_constant_scalar_space.hpp"


namespace Bempp
{

  %template(PiecewiseConstantScalarSpace_float64)     PiecewiseConstantScalarSpace<double>;
  %template(PiecewiseConstantScalarSpace_float32)     PiecewiseConstantScalarSpace<float>;
  %template(PiecewiseConstantScalarSpace_complex64)   PiecewiseConstantScalarSpace<std::complex<float> >;
  %template(PiecewiseConstantScalarSpace_complex128)  PiecewiseConstantScalarSpace<std::complex<double> >;

}


%pythoncode %{

def PiecewiseConstantScalarSpace(grid,BasisFunctionType='float64'):
    """Space of piecewise constant scalar functions"""
    name='PiecewiseConstantScalarSpace'+'_'+checkType(BasisFunctionType)
    return globals()[name](grid)
	 
 %}

