%{
#include "space/space.hpp"
#include <complex>
%}


// TODO
// %include "grid_docstrings.i"


%include "space/scalar_space.hpp";

namespace Bempp
{

  %template(ScalarSpace_float64)     ScalarSpace<double>;
  %template(ScalarSpace_float32)     ScalarSpace<float>;
  %template(ScalarSpace_complex64)   ScalarSpace<std::complex<float> >;
  %template(ScalarSpace_complex128)  ScalarSpace<std::complex<double> >;


} 

