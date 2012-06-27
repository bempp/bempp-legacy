%{
#include "assembly/assembly_options.hpp"
%}

namespace Bempp
{

    %extend AssemblyOptions
    {
        %ignore switchToFmm;
        %ignore switchToOpenCl;
        %ignore switchToTbb;
        %ignore parallelisationOptions;
    }

} // namespace Bempp

%include "assembly/assembly_options.hpp"
