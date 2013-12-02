%{
#include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "assembly/laplace_3d_calderon_projector.hpp"
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
    laplace3dHypersingularBoundaryOperatorPython;
%feature("compactdefaultargs")
    laplace3dExteriorCalderonProjector;
%feature("compactdefaultargs")
    laplace3dInteriorCalderonProjector;

%ignore laplace3dHypersingularBoundaryOperator;

} // namespace Bempp

%newobject laplace3dHypersingularBoundaryOperatorPython;

%inline %{

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dHypersingularBoundaryOperatorPython(
        const boost::shared_ptr<const Context<BasisFunctionType,ResultType> >& context,
        const boost::shared_ptr<const Space<BasisFunctionType> >& domain,
        const boost::shared_ptr<const Space<BasisFunctionType> >& range,
        const boost::shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY) {

    return laplace3dHypersingularBoundaryOperator(context,domain,range,dualToRange,
                                                  label,symmetry);


}

} // namespace Bempp


%}



#define shared_ptr boost::shared_ptr
%include "assembly/laplace_3d_single_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/laplace_3d_calderon_projector.hpp"
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
    laplace3dHypersingularBoundaryOperatorPython);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dExteriorCalderonProjector);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    laplace3dInteriorCalderonProjector);

} // namespace Bempp

