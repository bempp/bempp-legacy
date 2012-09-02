%{
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

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
        context, domain, range, dualToRange, waveNumber):
     """Construct a hypersingular operator for the modified Helmholtz equation in 3D."""
     return _constructModifiedHelmholtzOperator(
         "modifiedHelmholtz3dHypersingularBoundaryOperator", context, 
         domain, range, dualToRange, waveNumber)

%}
