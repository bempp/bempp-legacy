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

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dSingleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dDoubleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dAdjointDoubleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dHypersingularOperator);

} // namespace Bempp

%pythoncode %{

def laplace3dSingleLayerPotential(testSpace, trialSpace, resultType=None):
    """Construct a single-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dSingleLayerPotential", testSpace, trialSpace, resultType)

def laplace3dDoubleLayerPotential(testSpace, trialSpace, resultType=None):
    """Construct a double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dDoubleLayerPotential", testSpace, trialSpace, resultType)

def laplace3dAdjointDoubleLayerPotential(testSpace, trialSpace, resultType=None):
    """Construct an adjoint double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dAdjointDoubleLayerPotential", testSpace, trialSpace, resultType)

def laplace3dHypersingularOperator(testSpace, trialSpace, resultType=None):
    """Construct a hypersingular operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dHypersingularOperator", testSpace, trialSpace, resultType)

%}
