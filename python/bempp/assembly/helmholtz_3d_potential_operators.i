%{
#include "assembly/helmholtz_3d_potential_operator_base.hpp"
#include "assembly/helmholtz_3d_single_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_potential_operator.hpp"
%}


#define shared_ptr boost::shared_ptr
%include "assembly/helmholtz_3d_potential_operator_base.hpp"
%include "assembly/helmholtz_3d_single_layer_potential_operator.hpp"
%include "assembly/helmholtz_3d_double_layer_potential_operator.hpp"
#undef shared shared_ptr

%define BEMPP_INSTANTIATE_HELMHOLTZ_POTENTIAL_3D_BASE(BASIS, PY_BASIS)
    %template(Helmholtz3dPotentialOperatorBase_Single_ ## _ ## PY_BASIS)
        Helmholtz3dPotentialOperatorBase<
        Helmholtz3dSingleLayerPotentialOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dPotentialOperatorBase_Double_ ## _ ## PY_BASIS)
        Helmholtz3dPotentialOperatorBase<
        Helmholtz3dDoubleLayerPotentialOperatorImpl< BASIS >, BASIS >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dSingleLayerPotentialOperatorImpl);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dDoubleLayerPotentialOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_TYPES(BEMPP_INSTANTIATE_HELMHOLTZ_POTENTIAL_3D_BASE);

} // namespace Bempp

