%{
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs")
    modifiedHelmholtz3dSingleLayerBoundaryOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dHypersingularBoundaryOperator;
} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dHypersingularBoundaryOperator);

} // namespace Bempp
