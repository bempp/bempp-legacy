%{
#include "assembly/modified_helmholtz_3d_single_layer_potential_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_potential_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_potential_operator.hpp"
// #include "assembly/modified_helmholtz_3d_hypersingular_operator.hpp"
%}

// TODO
// %include "modified_helmholtz_3d_operators_docstrings.i"

namespace Bempp
{

%extend ModifiedHelmholtz3dSingleLayerPotentialOperator { %ignore clone; }
%extend ModifiedHelmholtz3dDoubleLayerPotentialOperator { %ignore clone; }
%extend ModifiedHelmholtz3dAdjointDoubleLayerPotentialOperator { %ignore clone; }
// %extend ModifiedHelmholtz3dHypersingularOperator { %ignore clone; }

} // namespace Bempp

%include "assembly/modified_helmholtz_3d_operator_base.hpp"
%include "assembly/modified_helmholtz_3d_single_layer_potential_operator.hpp"
%include "assembly/modified_helmholtz_3d_double_layer_potential_operator.hpp"
%include "assembly/modified_helmholtz_3d_adjoint_double_layer_potential_operator.hpp"
// %include "assembly/modified_helmholtz_3d_hypersingular_operator.hpp"

%define BEMPP_INSTANTIATE_MODIFIED_HELMHOLTZ_3D_BASE(BASIS, KERNEL, RESULT, PY_BASIS, PY_KERNEL, PY_RESULT)
    %template(ModifiedHelmholtz3dOperatorBase_Single_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dOperatorBase<
        ModifiedHelmholtz3dSingleLayerPotentialOperatorImpl< BASIS, KERNEL, RESULT >, BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dOperatorBase_Double_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dOperatorBase<
        ModifiedHelmholtz3dDoubleLayerPotentialOperatorImpl< BASIS, KERNEL, RESULT >, BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dOperatorBase_AdjointDouble_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dOperatorBase<
        ModifiedHelmholtz3dAdjointDoubleLayerPotentialOperatorImpl< BASIS, KERNEL, RESULT >, BASIS, KERNEL, RESULT >;

//    %template(ModifiedHelmholtz3dOperatorBase_Hypersingular_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
//        ModifiedHelmholtz3dOperatorBase<
//        ModifiedHelmholtz3dHypersingularOperatorImpl< BASIS, KERNEL, RESULT >, BASIS, KERNEL, RESULT >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dSingleLayerPotentialOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(BEMPP_INSTANTIATE_MODIFIED_HELMHOLTZ_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dDoubleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dAdjointDoubleLayerPotentialOperator);
// BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
//     ModifiedHelmholtz3dHypersingularOperator);

} // namespace Bempp

%pythoncode %{

def _constructModifiedHelmholtzOperator(className, domain, range, dualToRange, waveNumber, resultType):
    basisFunctionType = domain.basisFunctionType()
    if (basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of all spaces must be the same")
    waveNumberIsComplex = complex(waveNumber).imag != 0

    # determine resultType
    if waveNumberIsComplex and resultType in ("float32", "float64"):
        raise TypeError("Real result type given for a complex wave number")
    if resultType is None:
        if waveNumberIsComplex:
            resultType = "complex128"
        else:
            resultType = "float64"
    else:
        resultType = checkType(resultType)

    # determine kernelType
    if waveNumberIsComplex:
        kernelType = resultType
    else:
        if resultType in ("float32", "complex64"):
            kernelType = "float32"
        else:
            kernelType = "float64"

    # construct object
    result = constructObjectTemplatedOnBasisKernelAndResult(
        className, basisFunctionType, kernelType, resultType,
        domain, range, dualToRange, waveNumber)
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def modifiedHelmholtz3dSingleLayerPotentialOperator(
        domain, range, dualToRange, waveNumber, resultType=None):
    """Construct a single-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dSingleLayerPotentialOperator", domain, range, dualToRange,
        waveNumber, resultType)

def modifiedHelmholtz3dDoubleLayerPotentialOperator(
        domain, range, dualToRange, waveNumber, resultType=None):
    """Construct a double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dDoubleLayerPotentialOperator", domain, range, dualToRange,
        waveNumber, resultType)

def modifiedHelmholtz3dAdjointDoubleLayerPotentialOperator(
        domain, range, dualToRange, waveNumber, resultType=None):
    """Construct an adjoint double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dAdjointDoubleLayerPotentialOperator", domain, range, dualToRange,
        waveNumber, resultType)

# def modifiedHelmholtz3dHypersingularOperator(
#         domain, range, dualToRange, waveNumber, resultType=None):
#     """Construct a hypersingular operator for the modified Helmholtz equation in 3D."""
#     return _constructModifiedHelmholtzOperator(
#         "ModifiedHelmholtz3dHypersingularOperator", domain, range, dualToRange,
#         waveNumber, resultType)

%}
