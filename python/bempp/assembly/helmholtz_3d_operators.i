%{
#include "assembly/helmholtz_3d_boundary_operator_base.hpp"
#include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

// TODO
// %include "helmholtz_3d_operators_docstrings.i"

#define shared_ptr boost::shared_ptr
%include "assembly/helmholtz_3d_boundary_operator_base.hpp"
%include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
%include "assembly/helmholtz_3d_double_layer_boundary_operator.hpp"
%include "assembly/helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
#undef shared_ptr

%define BEMPP_INSTANTIATE_HELMHOLTZ_3D_BASE(BASIS, PY_BASIS)
    %template(Helmholtz3dBoundaryOperatorBase_Single_ ## _ ## PY_BASIS)
        Helmholtz3dBoundaryOperatorBase<
        Helmholtz3dSingleLayerBoundaryOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dBoundaryOperatorBase_Double_ ## _ ## PY_BASIS)
        Helmholtz3dBoundaryOperatorBase<
        Helmholtz3dDoubleLayerBoundaryOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dBoundaryOperatorBase_AdjointDouble_ ## _ ## PY_BASIS)
        Helmholtz3dBoundaryOperatorBase<
        Helmholtz3dAdjointDoubleLayerBoundaryOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dBoundaryOperatorBase_Hypersingular_ ## _ ## PY_BASIS)
        Helmholtz3dBoundaryOperatorBase<
        Helmholtz3dHypersingularBoundaryOperatorImpl< BASIS >, BASIS >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dSingleLayerBoundaryOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_TYPES(BEMPP_INSTANTIATE_HELMHOLTZ_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dHypersingularBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dHypersingularBoundaryOperator);

} // namespace Bempp

%pythoncode %{

def _constructHelmholtzOperator(className, context, domain, range, dualToRange, waveNumber):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of context and all spaces must be the same")
    resultType = context.resultType()
    result = constructObjectTemplatedOnBasis(
        className, basisFunctionType, context, domain, range, dualToRange, waveNumber)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def helmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def helmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def helmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def helmholtz3dHypersingularBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dHypersingularBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

%}
