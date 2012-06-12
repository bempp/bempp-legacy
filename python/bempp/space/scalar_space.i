%{
#include "space/space.hpp"
%}


// TODO
// %include "grid_docstrings.i"


%include "space/scalar_space.hpp";

namespace Bempp
{
  %template(ScalarSpaceDouble) ScalarSpace<double>;
} 
