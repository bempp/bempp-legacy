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

%template(PiecewiseConstantScalarSpaceDouble) PiecewiseConstantScalarSpace<double>;

}
