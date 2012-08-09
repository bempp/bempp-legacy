%{
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

// TODO
// %include "modified_helmholtz_3d_operators_docstrings.i"

#define shared_ptr boost::shared_ptr
%include "assembly/modified_helmholtz_3d_boundary_operator_base.hpp"
%include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#undef shared_ptr

%define BEMPP_INSTANTIATE_MODIFIED_HELMHOLTZ_3D_BASE(BASIS, KERNEL, RESULT, PY_BASIS, PY_KERNEL, PY_RESULT)
    %template(ModifiedHelmholtz3dBoundaryOperatorBase_Single_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dSingleLayerBoundaryOperatorImpl< BASIS , KERNEL , RESULT >,
        BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dBoundaryOperatorBase_Double_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl< BASIS , KERNEL , RESULT >,
        BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dBoundaryOperatorBase_AdjointDouble_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperatorImpl< BASIS , KERNEL , RESULT >,
        BASIS, KERNEL, RESULT >;

    %template(ModifiedHelmholtz3dBoundaryOperatorBase_Hypersingular_ ## _ ## PY_BASIS ## _ ## PY_KERNEL ## _ ## PY_RESULT)
        ModifiedHelmholtz3dBoundaryOperatorBase<
        ModifiedHelmholtz3dHypersingularBoundaryOperatorImpl< BASIS , KERNEL , RESULT >,
        BASIS , KERNEL , RESULT >;
%enddef

namespace Bempp
{

  BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
  ModifiedHelmholtz3dSingleLayerBoundaryOperatorImpl);
  BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
  ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl);
  BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
  ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperatorImpl);
  BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
  ModifiedHelmholtz3dHypersingularBoundaryOperatorImpl);

BEMPP_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(BEMPP_INSTANTIATE_MODIFIED_HELMHOLTZ_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ModifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator);
 BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
     ModifiedHelmholtz3dHypersingularBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dHypersingularBoundaryOperator);

} // namespace Bempp

%pythoncode %{

def _constructModifiedHelmholtzOperator(className, context,
                                        domain, range, dualToRange, waveNumber):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of context and all spaces must be the same")
    resultType = context.resultType()

    waveNumberIsComplex = complex(waveNumber).imag != 0
    if waveNumberIsComplex and resultType in ("float32", "float64"):
        raise TypeError("Real result type given for a complex wave number")

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
        context, domain, range, dualToRange, waveNumber)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def modifiedHelmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a single-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def modifiedHelmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct an adjoint double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def modifiedHelmholtz3dHypersingularBoundaryOperator(
         domain, range, dualToRange, waveNumber):
     """Construct a hypersingular operator for the modified Helmholtz equation in 3D."""
     return _constructModifiedHelmholtzOperator(
         "modifiedHelmholtz3dHypersingularBoundaryOperator", context, domain, range, dualToRange,
         waveNumber)

%}
