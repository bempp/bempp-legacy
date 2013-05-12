%{
#include "assembly/assembly_options.hpp"
%}

namespace Bempp
{

%extend AssemblyOptions
{
    %ignore switchToTbb;
    %feature("compactdefaultargs") enableSingularIntegralCaching;
    %feature("compactdefaultargs") enableSparseStorageOfMassMatrices;
}

} // namespace Bempp

%include "assembly/assembly_options.hpp"
