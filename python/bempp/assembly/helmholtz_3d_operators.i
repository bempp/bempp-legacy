%{
#include "assembly/helmholtz_3d_boundary_operator_base.hpp"
#include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_hypersingular_boundary_operator.hpp"
%}

// Redeclared because SWIG doesn't parse correctly ...::KernelType.
// So we replace it with the explicit
// typename ScalarTraits<BasisFunctionType>::ComplexType

#define shared_ptr boost::shared_ptr
namespace Bempp
{

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
helmholtz3dSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);
template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
helmholtz3dDoubleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);
template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
helmholtz3dAdjointDoubleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);
template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
helmholtz3dHypersingularBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

%feature("compactdefaultargs") helmholtz3dSingleLayerBoundaryOperator;

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dSingleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dAdjointDoubleLayerBoundaryOperator);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(
    helmholtz3dHypersingularBoundaryOperator);

} // namespace Bempp

%pythoncode %{

def _constructHelmholtzOperator(className, context, domain, range, dualToRange, waveNumber):
    basisFunctionType = context.basisFunctionType()
    if (basisFunctionType != domain.basisFunctionType() or
            basisFunctionType != range.basisFunctionType() or
            basisFunctionType != dualToRange.basisFunctionType()):
        raise TypeError("BasisFunctionType of context and all spaces must be the same")
    resultType = context.resultType()
    result = constructObjectTemplatedOnBasis(
        className, basisFunctionType, context, domain, range, dualToRange, waveNumber)
    result._context = context
    result._domain = domain
    result._range = range
    result._dualToRange = dualToRange
    return result

def helmholtz3dSingleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dSingleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def helmholtz3dDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def helmholtz3dAdjointDoubleLayerBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dAdjointDoubleLayerBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

def helmholtz3dHypersingularBoundaryOperator(
        context, domain, range, dualToRange, waveNumber):
    """Construct a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "helmholtz3dHypersingularBoundaryOperator", context,
        domain, range, dualToRange, waveNumber)

%}
#undef shared_ptr
