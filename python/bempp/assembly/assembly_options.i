%{
#include "assembly/assembly_options.hpp"
%}

namespace Bempp
{

%feature("autodoc", "eps -> float") AcaOptions::eps;
%feature("autodoc", "eta -> float") AcaOptions::eta;
%feature("autodoc", "globalAssemblyBeforeCompression -> bool") AcaOptions::globalAssemblyBeforeCompression;
%feature("autodoc", "maximumBlockSize -> int") AcaOptions::maximumBlockSize;
%feature("autodoc", "maximumRank -> int") AcaOptions::maximumRank;
%feature("autodoc", "minimumBlockSize -> int") AcaOptions::minimumBlockSize;
%feature("autodoc", "outputFname -> string") AcaOptions::outputFname;
%feature("autodoc", "outputPostscript -> bool") AcaOptions::outputPostscript;
%feature("autodoc", "recompress -> bool") AcaOptions::recompress;
%feature("autodoc", "scaling -> float") AcaOptions::scaling;

%extend AssemblyOptions
{
    %ignore switchToTbb;
    %feature("compactdefaultargs") enableSingularIntegralCaching;
    %feature("compactdefaultargs") enableSparseStorageOfMassMatrices;
}

} // namespace Bempp

%include "assembly/aca_options.hpp"
%include "assembly/assembly_options.hpp"
