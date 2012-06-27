%{
#include "assembly/evaluation_options.hpp"
%}

namespace Bempp
{

    %extend EvaluationOptions
    {
        %ignore switchToOpenCl;
    }

} // namespace Bempp

%include "assembly/evaluation_options.hpp"
