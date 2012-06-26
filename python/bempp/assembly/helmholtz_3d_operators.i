%{
#include "assembly/helmholtz_3d_operator_base.hpp"
#include "assembly/helmholtz_3d_single_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_adjoint_double_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_operator.hpp"
%}

// TODO
// %include "helmholtz_3d_operators_docstrings.i"

namespace Bempp
{

%extend Helmholtz3dSingleLayerPotentialOperator { %ignore clone; }
%extend Helmholtz3dDoubleLayerPotentialOperator { %ignore clone; }
%extend Helmholtz3dAdjointDoubleLayerPotentialOperator { %ignore clone; }
%extend Helmholtz3dHypersingularOperator { %ignore clone; }

} // namespace Bempp

%include "assembly/helmholtz_3d_operator_base.hpp"
%include "assembly/helmholtz_3d_single_layer_potential_operator.hpp"
%include "assembly/helmholtz_3d_double_layer_potential_operator.hpp"
%include "assembly/helmholtz_3d_adjoint_double_layer_potential_operator.hpp"
%include "assembly/helmholtz_3d_hypersingular_operator.hpp"

%define BEMPP_INSTANTIATE_HELMHOLTZ_3D_BASE(BASIS, PY_BASIS)
    %template(Helmholtz3dOperatorBase_Single_ ## _ ## PY_BASIS)
        Helmholtz3dOperatorBase<
        Helmholtz3dSingleLayerPotentialOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dOperatorBase_Double_ ## _ ## PY_BASIS)
        Helmholtz3dOperatorBase<
        Helmholtz3dDoubleLayerPotentialOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dOperatorBase_AdjointDouble_ ## _ ## PY_BASIS)
        Helmholtz3dOperatorBase<
        Helmholtz3dAdjointDoubleLayerPotentialOperatorImpl< BASIS >, BASIS >;

    %template(Helmholtz3dOperatorBase_Hypersingular_ ## _ ## PY_BASIS)
        Helmholtz3dOperatorBase<
        Helmholtz3dHypersingularOperatorImpl< BASIS >, BASIS >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dSingleLayerPotentialOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_TYPES(BEMPP_INSTANTIATE_HELMHOLTZ_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dDoubleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dAdjointDoubleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    Helmholtz3dHypersingularOperator);

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

def helmholtz3dSingleLayerPotentialOperator(domain, range, dualToRange, waveNumber):
    """Construct a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dSingleLayerPotentialOperator", domain, range, dualToRange, waveNumber)

def helmholtz3dDoubleLayerPotentialOperator(domain, range, dualToRange, waveNumber):
    """Construct a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dDoubleLayerPotentialOperator", domain, range, dualToRange, waveNumber)

def helmholtz3dAdjointDoubleLayerPotentialOperator(domain, range, dualToRange, waveNumber):
    """Construct an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dAdjointDoubleLayerPotentialOperator", domain, range, dualToRange, waveNumber)

def helmholtz3dHypersingularOperator(domain, range, dualToRange, waveNumber):
    """Construct a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dHypersingularOperator", domain, range, dualToRange, waveNumber)

%}
