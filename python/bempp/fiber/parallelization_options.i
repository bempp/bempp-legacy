%{
#include "fiber/parallelization_options.hpp"
  %}

namespace Fiber
{

  %extend ParallelizationOptions
  {
    %ignore switchToOpenCl;
    %ignore openClOptions;
  }

}

%include "fiber/parallelization_options.hpp"
