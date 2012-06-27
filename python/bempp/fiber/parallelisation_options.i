%{
#include "fiber/parallelisation_options.hpp"
%}

namespace Fiber
{

    %extend ParallelisationOptions
    {
        %ignore switchToOpenCl;
        %ignore openClOptions;
    }

}

%include "fiber/parallelisation_options.hpp"
