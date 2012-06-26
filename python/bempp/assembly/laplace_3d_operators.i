%{
#include "assembly/laplace_3d_operator_base.hpp"
#include "assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_double_layer_potential_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_potential_operator.hpp"
#include "assembly/laplace_3d_hypersingular_operator.hpp"
%}

// TODO
// %include "laplace_3d_operators_docstrings.i"

namespace Bempp
{

// for some reason, SWIG doesn't see that this is an override of LinearOperator::clone()
// (which it has been told to ignore)
%extend Laplace3dSingleLayerPotentialOperator { %ignore clone; }
%extend Laplace3dDoubleLayerPotentialOperator { %ignore clone; }
%extend Laplace3dAdjointDoubleLayerPotentialOperator { %ignore clone; }
%extend Laplace3dHypersingularOperator { %ignore clone; }

}

%include "assembly/laplace_3d_operator_base.hpp"
%include "assembly/laplace_3d_single_layer_potential_operator.hpp"
%include "assembly/laplace_3d_double_layer_potential_operator.hpp"
%include "assembly/laplace_3d_adjoint_double_layer_potential_operator.hpp"
%include "assembly/laplace_3d_hypersingular_operator.hpp"

%define BEMPP_INSTANTIATE_LAPLACE_3D_BASE(BASIS, RESULT, PY_BASIS, PY_RESULT)
    %template(Laplace3dOperatorBase_Single_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dOperatorBase<
        Laplace3dSingleLayerPotentialOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;

    %template(Laplace3dOperatorBase_Double_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dOperatorBase<
        Laplace3dDoubleLayerPotentialOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;

    %template(Laplace3dOperatorBase_AdjointDouble_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dOperatorBase<
        Laplace3dAdjointDoubleLayerPotentialOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;

    %template(Laplace3dOperatorBase_Hypersingular_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dOperatorBase<
        Laplace3dHypersingularOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dSingleLayerPotentialOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_RESULT_TYPES(BEMPP_INSTANTIATE_LAPLACE_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dDoubleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dAdjointDoubleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dHypersingularOperator);

} // namespace Bempp

%pythoncode %{

def laplace3dSingleLayerPotentialOperator(domain, range, dualToRange, resultType=None):
    """Construct a single-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dSingleLayerPotentialOperator", domain, range, dualToRange, resultType)

def laplace3dDoubleLayerPotentialOperator(domain, range, dualToRange, resultType=None):
    """Construct a double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dDoubleLayerPotentialOperator", domain, range, dualToRange, resultType)

def laplace3dAdjointDoubleLayerPotentialOperator(domain, range, dualToRange, resultType=None):
    """Construct an adjoint double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dAdjointDoubleLayerPotentialOperator", domain, range, dualToRange, resultType)

def laplace3dHypersingularOperator(domain, range, dualToRange, resultType=None):
    """Construct a hypersingular operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dHypersingularOperator", domain, range, dualToRange, resultType)

%}
