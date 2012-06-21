%{
#include "assembly/surface_normal_dependent_functor.hpp"
%}

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(SurfaceNormalDependentFunctor);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(SurfaceNormalDependentFunctor);
} // namespace Bempp

%include "assembly/surface_normal_dependent_functor.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(SurfaceNormalDependentFunctor);
} // namespace Bempp
