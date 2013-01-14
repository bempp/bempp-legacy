%{
#include "assembly/maxwell_3d_single_layer_boundary_operator.hpp"
#include "assembly/maxwell_3d_double_layer_boundary_operator.hpp"
#include "assembly/maxwell_3d_identity_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs")
    maxwell3dSingleLayerBoundaryOperator;
%feature("compactdefaultargs")
    maxwell3dDoubleLayerBoundaryOperator;
%feature("compactdefaultargs")
    maxwell3dIdentityOperator;
} // namespace Bempp

// Redeclared because SWIG doesn't parse correctly ...::KernelType.
// So we replace it with the explicit
// typename ScalarTraits<BasisFunctionType>::ComplexType

#define shared_ptr boost::shared_ptr

%include "assembly/maxwell_3d_identity_operator.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY,
        bool useInterpolation = false,
        int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);
template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dDoubleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY,
        bool useInterpolation = false,
        int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    maxwell3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    maxwell3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
    maxwell3dIdentityOperator);

} // namespace Bempp
#undef shared_ptr
