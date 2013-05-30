%{
#include "assembly/modified_helmholtz_3d_single_layer_potential_operator.hpp"
#include "assembly/modified_helmholtz_3d_double_layer_potential_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs")
    modifiedHelmholtz3dSingleLayerPotentialOperator;
%feature("compactdefaultargs")
    modifiedHelmholtz3dDoubleLayerPotentialOperator;
} // namespace Bempp

// Swig isn't able to parse the
// Helmholtz3d...PotentialOperator<BasisFunctionType>::KernelType
// type -- replace it with the explicit
// ScalarTraits<BasisFunctionType>::ComplexType.

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
    boost::shared_ptr<PotentialOperator<
        BasisFunctionType,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
    >
modifiedHelmholtz3dSingleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::ModifiedHelmholtz3dSingleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

template <typename BasisFunctionType>
    boost::shared_ptr<PotentialOperator<
        BasisFunctionType,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
    >
modifiedHelmholtz3dDoubleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::ModifiedHelmholtz3dDoubleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(modifiedHelmholtz3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(modifiedHelmholtz3dDoubleLayerPotentialOperator);
}

