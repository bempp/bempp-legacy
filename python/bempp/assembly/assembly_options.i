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
    %feature("compactdefaultargs") enableJointAssembly;
}

} // namespace Bempp

%include "assembly/assembly_options.hpp"
