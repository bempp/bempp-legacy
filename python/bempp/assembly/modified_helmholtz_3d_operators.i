%{
#include "assembly/modified_helmholtz_3d_single_layer_potential.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_potential.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_potential.hpp"
// #include "assembly/modified_helmholtz_3d_hypersingular_operator.hpp"
%}

// TODO
// %include "modified_helmholtz_3d_operators_docstrings.i"

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dSingleLayerPotential);
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dDoubleLayerPotential);
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dAdjointDoubleLayerPotential);
// BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dHypersingularOperator);

} // namespace Bempp

%include "assembly/modified_helmholtz_3d_single_layer_potential.hpp"
%include "assembly/modified_helmholtz_3d_double_layer_potential.hpp"
%include "assembly/modified_helmholtz_3d_adjoint_double_layer_potential.hpp"
// %include "assembly/modified_helmholtz_3d_hypersingular_operator.hpp"

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dSingleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dDoubleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dAdjointDoubleLayerPotential);
// BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ModifiedHelmholtz3dHypersingularOperator);

} // namespace Bempp

%pythoncode %{

def _constructModifiedHelmholtzOperator(className, testSpace, trialSpace, waveNumber, resultType):
    basisFunctionType = testSpace.basisFunctionType()
    if (basisFunctionType != trialSpace.basisFunctionType()):
        raise TypeError("BasisFunctionType of testSpace must match that of trialSpace")
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
        testSpace, trialSpace, waveNumber)
    result._testSpace = testSpace
    result._trialSpace = trialSpace
    return result

def modifiedHelmholtz3dSingleLayerPotential(testSpace, trialSpace, waveNumber, resultType=None):
    """Construct a single-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dSingleLayerPotential", testSpace, trialSpace, waveNumber, resultType)

def modifiedHelmholtz3dDoubleLayerPotential(testSpace, trialSpace, waveNumber):
    """Construct a double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dDoubleLayerPotential", testSpace, trialSpace, waveNumber, resultType)

def modifiedHelmholtz3dAdjointDoubleLayerPotential(testSpace, trialSpace, waveNumber):
    """Construct an adjoint double-layer-potential operator for the modified Helmholtz equation in 3D."""
    return _constructModifiedHelmholtzOperator(
        "ModifiedHelmholtz3dAdjointDoubleLayerPotential", testSpace, trialSpace, waveNumber, resultType)

# def modifiedHelmholtz3dHypersingularOperator(testSpace, trialSpace, waveNumber):
#     """Construct a hypersingular operator for the modified Helmholtz equation in 3D."""
#     return _constructModifiedHelmholtzOperator(
#         "ModifiedHelmholtz3dHypersingularOperator", testSpace, trialSpace, waveNumber, resultType)

%}
