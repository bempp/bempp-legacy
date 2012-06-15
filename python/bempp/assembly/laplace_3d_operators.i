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
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dSingleLayerPotential);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dDoubleLayerPotential);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dAdjointDoubleLayerPotential);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dHypersingularOperator);
}

%pythoncode %{

class Laplace3dSingleLayerPotential(Template2, ElementarySingularIntegralOperator):
    """Single-layer-potential operator for the Laplace equation in 3D."""
    def __init__(self, basisFunctionType, resultType, *args, **kwargs):
        """Initialise operator."""
        super(Laplace3dSingleLayerPotential,self).\
            __init__('Laplace3dSingleLayerPotential', basisFunctionType, resultType,
                     *args, **kwargs)

class Laplace3dDoubleLayerPotential(Template2, ElementarySingularIntegralOperator):
    """Single-layer-potential operator for the Laplace equation in 3D."""
    def __init__(self, basisFunctionType, resultType, *args, **kwargs):
        """Initialise operator."""
        super(Laplace3dDoubleLayerPotential,self).\
            __init__('Laplace3dDoubleLayerPotential', basisFunctionType, resultType,
                     *args, **kwargs)

class Laplace3dAdjointDoubleLayerPotential(Template2, ElementarySingularIntegralOperator):
    """Single-layer-potential operator for the Laplace equation in 3D."""
    def __init__(self, basisFunctionType, resultType, *args, **kwargs):
        """Initialise operator."""
        super(Laplace3dAdjointDoubleLayerPotential,self).\
            __init__('Laplace3dAdjointDoubleLayerPotential', basisFunctionType, resultType,
                     *args, **kwargs)

class Laplace3dHypersingularOperator(Template2, ElementarySingularIntegralOperator):
    """Single-layer-potential operator for the Laplace equation in 3D."""
    def __init__(self, basisFunctionType, resultType, *args, **kwargs):
        """Initialise operator."""
        super(Laplace3dHypersingularOperator,self).\
            __init__('Laplace3dHypersingularOperator', basisFunctionType, resultType,
                     *args, **kwargs)

# TODO: docs

%}
