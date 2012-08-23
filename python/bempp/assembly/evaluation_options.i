%{
#include "assembly/evaluation_options.hpp"
%}

namespace Bempp
{

%extend EvaluationOptions
{
    %ignore switchToTbb;
}

} // namespace Bempp

%include "assembly/evaluation_options.hpp"

    
