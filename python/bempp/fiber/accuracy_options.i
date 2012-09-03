%{
#include "fiber/accuracy_options.hpp"
%}

namespace Fiber
{

%feature("autodoc", "doubleRegular -> QuadratureOptions") AccuracyOptions::doubleRegular;
%feature("autodoc", "doubleSingular -> QuadratureOptions") AccuracyOptions::doubleSingular;
%feature("autodoc", "singleRegular -> QuadratureOptions") AccuracyOptions::singleRegular;

} // namespace Fiber

%include "fiber/accuracy_options.hpp"
