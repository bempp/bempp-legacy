%{
#include "assembly/maxwell_3d_single_layer_potential_operator.hpp"
#include "assembly/maxwell_3d_double_layer_potential_operator.hpp"
// #include "assembly/maxwell_3d_far_field_single_layer_potential_operator.hpp"
// #include "assembly/maxwell_3d_far_field_double_layer_potential_operator.hpp"
%}

// Swig isn't able to parse the
// Maxwell3d...PotentialOperator<BasisFunctionType>::KernelType
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
maxwell3dSingleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::Maxwell3dSingleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

template <typename BasisFunctionType>
    boost::shared_ptr<PotentialOperator<
        BasisFunctionType,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
    >
maxwell3dDoubleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::Maxwell3dDoubleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

// template <typename BasisFunctionType>
//     boost::shared_ptr<PotentialOperator<
//         BasisFunctionType,
//         typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
//     >
// maxwell3dFarFieldSingleLayerPotentialOperator(
//     typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
// {
//     typedef Bempp::Maxwell3dFarFieldSingleLayerPotentialOperator<BasisFunctionType> Type;
//     return boost::shared_ptr<Type>(new Type(waveNumber));
// }

// template <typename BasisFunctionType>
//     boost::shared_ptr<PotentialOperator<
//         BasisFunctionType,
//         typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
//     >
// maxwell3dFarFieldDoubleLayerPotentialOperator(
//     typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
// {
//     typedef Bempp::Maxwell3dFarFieldDoubleLayerPotentialOperator<BasisFunctionType> Type;
//     return boost::shared_ptr<Type>(new Type(waveNumber));
// }

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(maxwell3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(maxwell3dDoubleLayerPotentialOperator);
// BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(maxwell3dFarFieldSingleLayerPotentialOperator);
// BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(maxwell3dFarFieldDoubleLayerPotentialOperator);
}

