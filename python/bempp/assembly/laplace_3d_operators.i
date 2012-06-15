%{
#include "assembly/laplace_3d_single_layer_potential.hpp"
#include "assembly/laplace_3d_double_layer_potential.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_potential.hpp"
#include "assembly/laplace_3d_hypersingular_operator.hpp"
%}

// TODO
// %include "laplace_3d_operators_docstrings.i"

%include "assembly/laplace_3d_single_layer_potential.hpp"
%include "assembly/laplace_3d_double_layer_potential.hpp"
%include "assembly/laplace_3d_adjoint_double_layer_potential.hpp"
%include "assembly/laplace_3d_hypersingular_operator.hpp"

namespace Bempp
{

BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dSingleLayerPotential);
BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dDoubleLayerPotential);
BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dAdjointDoubleLayerPotential);
BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dHypersingularOperator);

} // namespace Bempp

%pythoncode %{

def _constructLaplaceOperator(className, resultType, testSpace, trialSpace):
    basisFunctionType = testSpace._basisFunctionType
    if (basisFunctionType != trialSpace._basisFunctionType):
        raise TypeError("BasisFunctionType of testSpace must match that of trialSpace")
    return constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType, 
        testSpace, trialSpace)

def laplace3dSingleLayerPotential(resultType, testSpace, trialSpace):
    """Construct a single-layer-potential operator for the Laplace equation in 3D."""
    return _constructLaplaceOperator(
        "Laplace3dSingleLayerPotential", resultType, testSpace, trialSpace)

def laplace3dDoubleLayerPotential(resultType, testSpace, trialSpace):
    """Construct a double-layer-potential operator for the Laplace equation in 3D."""
    return _constructLaplaceOperator(
        "Laplace3dDoubleLayerPotential", resultType, testSpace, trialSpace)

def laplace3dAdjointDoubleLayerPotential(resultType, testSpace, trialSpace):
    """Construct an adjoint double-layer-potential operator for the Laplace equation in 3D."""
    return _constructLaplaceOperator(
        "Laplace3dAdjointDoubleLayerPotential", resultType, testSpace, trialSpace)

def laplace3dHypersingularOperator(resultType, testSpace, trialSpace):
    """Construct a hypersingular operator for the Laplace equation in 3D."""
    return _constructLaplaceOperator(
        "Laplace3dHypersingularOperator", resultType, testSpace, trialSpace)

%}
