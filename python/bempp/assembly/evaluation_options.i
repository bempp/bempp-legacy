%{
#include "assembly/evaluation_options.hpp"
  %}

namespace Bempp {

  %extend EvaluationOptions {
    
    %ignore switchToOpenCl;
    %feature("compactdefaultargs") switchToTbb;

  }

}

%include "assembly/evaluation_options.hpp"

    
