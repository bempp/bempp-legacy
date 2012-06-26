%{
#include "assembly/laplace_3d_boundary_operator_base.hpp"
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
%}

// TODO
// %include "laplace_3d_operators_docstrings.i"

namespace Bempp
{

// for some reason, SWIG doesn't see that this is an override of BoundaryOperator::clone()
// (which it has been told to ignore)
%extend Laplace3dSingleLayerBoundaryOperator { %ignore clone; }
%extend Laplace3dDoubleLayerBoundaryOperator { %ignore clone; }
%extend Laplace3dAdjointDoubleLayerBoundaryOperator { %ignore clone; }
%extend Laplace3dHypersingularBoundaryOperator { %ignore clone; }

}

%include "assembly/laplace_3d_boundary_operator_base.hpp"
%include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"

%define BEMPP_INSTANTIATE_LAPLACE_3D_BASE(BASIS, RESULT, PY_BASIS, PY_RESULT)
    %template(Laplace3dBoundaryOperatorBase_Single_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dBoundaryOperatorBase<
        Laplace3dSingleLayerBoundaryOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;

    %template(Laplace3dBoundaryOperatorBase_Double_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dBoundaryOperatorBase<
        Laplace3dDoubleLayerBoundaryOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;

    %template(Laplace3dBoundaryOperatorBase_AdjointDouble_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dBoundaryOperatorBase<
        Laplace3dAdjointDoubleLayerBoundaryOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;

    %template(Laplace3dBoundaryOperatorBase_Hypersingular_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dBoundaryOperatorBase<
        Laplace3dHypersingularBoundaryOperatorImpl< BASIS, RESULT >, BASIS, RESULT >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dSingleLayerBoundaryOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_AND_RESULT_TYPES(BEMPP_INSTANTIATE_LAPLACE_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dHypersingularBoundaryOperator);

} // namespace Bempp

%pythoncode %{

def laplace3dSingleLayerBoundaryOperator(domain, range, dualToRange, resultType=None):
    """Construct a single-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dSingleLayerBoundaryOperator", domain, range, dualToRange, resultType)

def laplace3dDoubleLayerBoundaryOperator(domain, range, dualToRange, resultType=None):
    """Construct a double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dDoubleLayerBoundaryOperator", domain, range, dualToRange, resultType)

def laplace3dAdjointDoubleLayerBoundaryOperator(domain, range, dualToRange, resultType=None):
    """Construct an adjoint double-layer-potential operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dAdjointDoubleLayerBoundaryOperator", domain, range, dualToRange, resultType)

def laplace3dHypersingularBoundaryOperator(domain, range, dualToRange, resultType=None):
    """Construct a hypersingular operator for the Laplace equation in 3D."""
    return _constructOperator(
    "Laplace3dHypersingularBoundaryOperator", domain, range, dualToRange, resultType)

%}
