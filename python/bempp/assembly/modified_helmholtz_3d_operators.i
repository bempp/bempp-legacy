%{
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_calderon_projector.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs")
    modifiedHelmholtz3dSingleLayerBoundaryOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dHypersingularBoundaryOperatorPython;
%feature("compactdefaultargs")
    modifiedHelmholtz3dExteriorCalderonProjector;


%ignore modifiedHelmholtz3dHypersingularBoundaryOperator;

} // namespace Bempp

%newobject modifiedHelmholtz3dHypersingularBoundaryOperatorPython;

%inline %{

namespace Bempp {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dHypersingularBoundaryOperatorPython(
        const boost::shared_ptr<const Context<BasisFunctionType,ResultType> >& context,
        const boost::shared_ptr<const Space<BasisFunctionType> >& domain,
        const boost::shared_ptr<const Space<BasisFunctionType> >& range,
        const boost::shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        KernelType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY,
        bool useInterpolation = false,
        int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY) {

    return modifiedHelmholtz3dHypersingularBoundaryOperator(context,domain,range,dualToRange,
                                                     waveNumber,label,symmetry,
                                                     useInterpolation,interpPtsPerWavelength);

}

} // namespace Bempp


%}


#define shared_ptr boost::shared_ptr
%include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
%include "assembly/modified_helmholtz_3d_calderon_projector.hpp"
//%include "assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
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
    modifiedHelmholtz3dHypersingularBoundaryOperatorPython);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    modifiedHelmholtz3dExteriorCalderonProjector);


} // namespace Bempp

