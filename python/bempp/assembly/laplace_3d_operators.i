%{
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs")
    laplace3dSingleLayerBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dAdjointDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dHypersingularBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dSyntheticSingleLayerBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dSyntheticDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dSyntheticAdjointDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    laplace3dSyntheticHypersingularBoundaryOperator;

// Do not emit warnings about ignored assignment operators.
%warnfilter(362) Laplace3dBoundaryOperatorBase::operator=;

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#undef shared_ptr

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dHypersingularBoundaryOperator);

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dSyntheticSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dSyntheticDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dSyntheticAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dSyntheticHypersingularBoundaryOperator);

} // namespace Bempp

