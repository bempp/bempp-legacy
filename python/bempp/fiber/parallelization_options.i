%{
#include "fiber/parallelization_options.hpp"
%}

namespace Fiber
{

%extend ParallelizationOptions
{
%feature("compactdefaultargs") setMaxThreadCount;
}

}

%include "fiber/parallelization_options.hpp"
