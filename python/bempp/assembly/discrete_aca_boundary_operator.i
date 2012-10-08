#ifdef WITH_AHMED

%{
#include "assembly/discrete_aca_boundary_operator.hpp"
#include "fiber/scalar_traits.hpp"
#include "common/shared_ptr.hpp"
%}

namespace Bempp 
{
    %ignore DiscreteAcaBoundaryOperator;
}

#define shared_ptr boost::shared_ptr
%include "assembly/discrete_aca_boundary_operator.hpp"
#undef boost::shared_ptr

%define BEMPP_ACA_FREE_FUNCTIONS(VALUE,PY_VALUE)
namespace Bempp 
{
    %template(acaOperatorApproximateLuInverse_## PY_VALUE) 
        acaOperatorApproximateLuInverse< VALUE >;
    %template(scaledAcaOperator_## PY_VALUE)
         scaledAcaOperator< VALUE >;
    %template(acaOperatorSum_## PY_VALUE) 
         acaOperatorSum< VALUE >;
}
%enddef

BEMPP_ITERATE_OVER_BASIS_TYPES(BEMPP_ACA_FREE_FUNCTIONS)

#endif // WITH_AHMED
