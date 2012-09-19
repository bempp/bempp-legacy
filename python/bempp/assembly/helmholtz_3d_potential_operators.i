%{
#include "assembly/helmholtz_3d_single_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_far_field_single_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_far_field_double_layer_potential_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs")
    helmholtz3dSingleLayerPotentialOperator;
%feature("compactdefaultargs")
    helmholtz3dDoubleLayerPotentialOperator;
%feature("compactdefaultargs")
    helmholtz3dFarFieldSingleLayerPotentialOperator;
%feature("compactdefaultargs")
    helmholtz3dFarFieldDoubleLayerPotentialOperator;
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
helmholtz3dSingleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::Helmholtz3dSingleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

template <typename BasisFunctionType>
    boost::shared_ptr<PotentialOperator<
        BasisFunctionType,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
    >
helmholtz3dDoubleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::Helmholtz3dDoubleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

template <typename BasisFunctionType>
    boost::shared_ptr<PotentialOperator<
        BasisFunctionType,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
    >
helmholtz3dFarFieldSingleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::Helmholtz3dFarFieldSingleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

template <typename BasisFunctionType>
    boost::shared_ptr<PotentialOperator<
        BasisFunctionType,
        typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType>
    >
helmholtz3dFarFieldDoubleLayerPotentialOperator(
    typename Bempp::ScalarTraits<BasisFunctionType>::ComplexType waveNumber)
{
    typedef Bempp::Helmholtz3dFarFieldDoubleLayerPotentialOperator<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(waveNumber));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(helmholtz3dSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(helmholtz3dDoubleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(helmholtz3dFarFieldSingleLayerPotentialOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(helmholtz3dFarFieldDoubleLayerPotentialOperator);
}

