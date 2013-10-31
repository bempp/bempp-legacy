%include "std_vector.i"

%{
#include "fiber/accuracy_options.hpp"
%}

namespace Fiber
{

%feature("autodoc", "doubleRegular -> QuadratureOptions") AccuracyOptions::doubleRegular;
%feature("autodoc", "doubleSingular -> QuadratureOptions") AccuracyOptions::doubleSingular;
%feature("autodoc", "singleRegular -> QuadratureOptions") AccuracyOptions::singleRegular;

%extend AccuracyOptionsEx
{
    %feature("compactdefaultargs") setDoubleRegular;
    %feature("compactdefaultargs") setDoubleSingular;
    %feature("compactdefaultargs") setSingleRegular;
}

} // namespace Fiber

%include "fiber/accuracy_options.hpp"
