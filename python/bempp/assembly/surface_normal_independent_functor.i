%{
#include "assembly/surface_normal_independent_functor.hpp"
%}

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(SurfaceNormalIndependentFunctor);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(SurfaceNormalIndependentFunctor);
} // namespace Bempp

%include "assembly/surface_normal_independent_functor.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(SurfaceNormalIndependentFunctor);
} // namespace Bempp
