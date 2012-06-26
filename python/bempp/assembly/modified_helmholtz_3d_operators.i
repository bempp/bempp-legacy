%{
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
// #include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

// TODO
// %include "modified_helmholtz_3d_operators_docstrings.i"

namespace Bempp
{

%extend ModifiedHelmholtz3dSingleLayerBoundaryOperator { %ignore clone; }
%extend ModifiedHelmholtz3dDoubleLayerBoundaryOperator { %ignore clone; }
%extend ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator { %ignore clone; }
// %extend ModifiedHelmholtz3dHypersingularBoundaryOperator { %ignore clone; }

} // namespace Bempp

%include "assembly/modified_helmholtz_3d_boundary_operator_base.hpp"
%include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
// %include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"

%define BEMPP_INSTANTIATE_MODIFIED_HELMHOLTZ_3D_BASE(BASIS, KERNEL, RESULT, PY_BASIS, PY_KERNEL, PY_RESULT)
    %template(ModifiedHelmholtz3dBoundaryOperatorBase_Single_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dSingleLayerBoundaryOperatorImpl< BASIS, KERNEL, RESULT >,
        BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dBoundaryOperatorBase_Double_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl< BASIS, KERNEL, RESULT >,
        BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dBoundaryOperatorBase_AdjointDouble_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperatorImpl< BASIS, KERNEL, RESULT >,
        BASIS, KERNEL, RESULT >;

//    %template(ModifiedHelmholtz3dBoundaryOperatorBase_Hypersingular_
//        ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
//        ModifiedHelmholtz3dBoundaryOperatorBase<
//        ModifiedHelmholtz3dHypersingularBoundaryOperatorImpl< BASIS, KERNEL, RESULT >,
//        BASIS, KERNEL, RESULT >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dSingleLayerBoundaryOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(BEMPP_INSTANTIATE_MODIFIED_HELMHOLTZ_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator);
// BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
//     ModifiedHelmholtz3dHypersingularBoundaryOperator);

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

def modifiedHelmholtz3dSingleLayerBoundaryOperator(
        domain, range, dualToRange, waveNumber, resultType=None):
    """Construct a single-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dSingleLayerBoundaryOperator", domain, range, dualToRange,
        waveNumber, resultType)

def modifiedHelmholtz3dDoubleLayerBoundaryOperator(
        domain, range, dualToRange, waveNumber, resultType=None):
    """Construct a double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dDoubleLayerBoundaryOperator", domain, range, dualToRange,
        waveNumber, resultType)

def modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator(
        domain, range, dualToRange, waveNumber, resultType=None):
    """Construct an adjoint double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator", domain, range, dualToRange,
        waveNumber, resultType)

# def modifiedHelmholtz3dHypersingularBoundaryOperator(
#         domain, range, dualToRange, waveNumber, resultType=None):
#     """Construct a hypersingular operator for the modified Helmholtz equation in 3D."""
#     return _constructModifiedHelmholtzOperator(
#         "ModifiedHelmholtz3dHypersingularBoundaryOperator", domain, range, dualToRange,
#         waveNumber, resultType)

%}
