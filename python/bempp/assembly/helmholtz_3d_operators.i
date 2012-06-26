%{
#include "assembly/helmholtz_3d_boundary_operator_base.hpp"
#include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

// TODO
// %include "helmholtz_3d_operators_docstrings.i"

namespace Bempp
{

%extend Helmholtz3dSingleLayerBoundaryOperator { %ignore clone; }
%extend Helmholtz3dDoubleLayerBoundaryOperator { %ignore clone; }
%extend Helmholtz3dAdjointDoubleLayerBoundaryOperator { %ignore clone; }
%extend Helmholtz3dHypersingularBoundaryOperator { %ignore clone; }

} // namespace Bempp

%include "assembly/helmholtz_3d_boundary_operator_base.hpp"
%include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
%include "assembly/helmholtz_3d_double_layer_boundary_operator.hpp"
%include "assembly/helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"

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

} // namespace Bempp

%pythoncode %{

def _constructHelmholtzOperator(className, domain, range, dualToRange, waveNumber):
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")
    resultType = promoteTypeToComplex(basisFunctionType)
    result = constructObjectTemplatedOnBasis(
        className, basisFunctionType, domain, range, dualToRange, waveNumber)
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def helmholtz3dSingleLayerBoundaryOperator(domain, range, dualToRange, waveNumber):
    """Construct a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dSingleLayerBoundaryOperator", domain, range, dualToRange, waveNumber)

def helmholtz3dDoubleLayerBoundaryOperator(domain, range, dualToRange, waveNumber):
    """Construct a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dDoubleLayerBoundaryOperator", domain, range, dualToRange, waveNumber)

def helmholtz3dAdjointDoubleLayerBoundaryOperator(domain, range, dualToRange, waveNumber):
    """Construct an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dAdjointDoubleLayerBoundaryOperator", domain, range, dualToRange, waveNumber)

def helmholtz3dHypersingularBoundaryOperator(domain, range, dualToRange, waveNumber):
    """Construct a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dHypersingularBoundaryOperator", domain, range, dualToRange, waveNumber)

%}
