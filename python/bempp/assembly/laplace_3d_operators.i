%{
#include "assembly/laplace_3d_boundary_operator_base.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
%}

#define shared_ptr boost::shared_ptr
%include "assembly/laplace_3d_boundary_operator_base.hpp"
%include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dHypersingularBoundaryOperator);

} // namespace Bempp

%pythoncode %{

def laplace3dSingleLayerBoundaryOperator(context, domain, range, dualToRange):
    """Construct a single-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dSingleLayerBoundaryOperator", context, domain, range, dualToRange)

def laplace3dDoubleLayerBoundaryOperator(context, domain, range, dualToRange):
    """Construct a double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dDoubleLayerBoundaryOperator", context, domain, range, dualToRange)

def laplace3dAdjointDoubleLayerBoundaryOperator(context, domain, range, dualToRange):
    """Construct an adjoint double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dAdjointDoubleLayerBoundaryOperator", context, domain, range, dualToRange)

def laplace3dHypersingularBoundaryOperator(context, domain, range, dualToRange):
    """Construct a hypersingular operator for the Laplace equation in 3D."""
    return _constructOperator(
    "laplace3dHypersingularBoundaryOperator", context, domain, range, dualToRange)

%}

