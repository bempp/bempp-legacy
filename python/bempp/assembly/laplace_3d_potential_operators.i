%{
#include "assembly/laplace_3d_potential_operator_base.hpp"
#include "assembly/laplace_3d_single_layer_potential_operator.hpp"
#include "assembly/laplace_3d_double_layer_potential_operator.hpp"
%}

namespace Bempp {


  BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Laplace3dSingleLayerPotentialOperator)
  BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Laplace3dDoubleLayerPotentialOperator)

  BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS(Laplace3dSingleLayerPotentialOperator)
  BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS(Laplace3dDoubleLayerPotentialOperator)
    }


#define shared_ptr boost::shared_ptr
%include "assembly/laplace_3d_potential_operator_base.hpp"
%include "assembly/laplace_3d_single_layer_potential_operator.hpp"
%include "assembly/laplace_3d_double_layer_potential_operator.hpp"
#undef shared_ptr

%define BEMPP_INSTANTIATE_LAPLACE_POTENTIAL_3D_BASE(BASIS, RESULT,  PY_BASIS, PY_RESULT)
%template(Laplace3dPotentialOperatorBase_Single_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dPotentialOperatorBase<
Laplace3dSingleLayerPotentialOperatorImpl< BASIS , RESULT >, BASIS , RESULT >;

    %template(Laplace3dPotentialOperatorBase_Double_ ## _ ## PY_BASIS ## _ ## PY_RESULT)
        Laplace3dPotentialOperatorBase<
    Laplace3dDoubleLayerPotentialOperatorImpl< BASIS , RESULT >, BASIS , RESULT >;
%enddef

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dSingleLayerPotentialOperatorImpl);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dDoubleLayerPotentialOperatorImpl);
BEMPP_ITERATE_OVER_BASIS_AND_RESULT_TYPES(BEMPP_INSTANTIATE_LAPLACE_POTENTIAL_3D_BASE);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    Laplace3dDoubleLayerPotentialOperator);



} // namespace Bempp

%pythoncode %{

def _constructLaplacePotentialOperator(className, context):
    basisFunctionType = context.basisFunctionType()
    resultType = context.resultType()
    result = constructObjectTemplatedOnBasisAndResult(
        className, basisFunctionType, resultType )
    result._context = context
    return result

	  %}
