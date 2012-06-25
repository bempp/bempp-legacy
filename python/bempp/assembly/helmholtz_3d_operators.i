%{
#include "assembly/helmholtz_3d_single_layer_potential.hpp"
#include "assembly/helmholtz_3d_double_layer_potential.hpp"
#include "assembly/helmholtz_3d_adjoint_double_layer_potential.hpp"
#include "assembly/helmholtz_3d_hypersingular_operator.hpp"
%}

// TODO
// %include "helmholtz_3d_operators_docstrings.i"

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dSingleLayerPotential);
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dDoubleLayerPotential);
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotential);
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);

// Workaround for SWIG being unable to decipher
// "typename ScalarTraits<...>::ComplexType"

%extend Helmholtz3dSingleLayerPotential
{
    Helmholtz3dSingleLayerPotential(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::complex<double>& waveNumber) {
            return new Bempp::Helmholtz3dSingleLayerPotential<BasisFunctionType >(
                testSpace, trialSpace,
                static_cast<Bempp::Helmholtz3dSingleLayerPotential<BasisFunctionType >::KernelType>(waveNumber));
        }
}

%extend Helmholtz3dSingleLayerPotential
{
    %ignore Helmholtz3dSingleLayerPotential;
}

%extend Helmholtz3dDoubleLayerPotential
{
    Helmholtz3dDoubleLayerPotential(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::complex<double>& waveNumber) {
            return new Bempp::Helmholtz3dDoubleLayerPotential<BasisFunctionType >(
                testSpace, trialSpace,
                static_cast<Bempp::Helmholtz3dDoubleLayerPotential<BasisFunctionType >::KernelType>(waveNumber));
        }
}

%extend Helmholtz3dDoubleLayerPotential
{
    %ignore Helmholtz3dDoubleLayerPotential;
}

%extend Helmholtz3dAdjointDoubleLayerPotential
{
    Helmholtz3dAdjointDoubleLayerPotential(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::complex<double>& waveNumber) {
            return new Bempp::Helmholtz3dAdjointDoubleLayerPotential<BasisFunctionType >(
                testSpace, trialSpace,
                static_cast<Bempp::Helmholtz3dAdjointDoubleLayerPotential<BasisFunctionType >::KernelType>(waveNumber));
        }
}

%extend Helmholtz3dAdjointDoubleLayerPotential
{
    %ignore Helmholtz3dAdjointDoubleLayerPotential;
}

%extend Helmholtz3dHypersingularOperator
{
    Helmholtz3dHypersingularOperator(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::complex<double>& waveNumber) {
            return new Bempp::Helmholtz3dHypersingularOperator<BasisFunctionType >(
                testSpace, trialSpace,
                static_cast<Bempp::Helmholtz3dHypersingularOperator<BasisFunctionType >::KernelType>(waveNumber));
        }
}

%extend Helmholtz3dHypersingularOperator
{
    %ignore Helmholtz3dHypersingularOperator;
}

} // namespace Bempp

%include "assembly/helmholtz_3d_single_layer_potential.hpp"
%include "assembly/helmholtz_3d_double_layer_potential.hpp"
%include "assembly/helmholtz_3d_adjoint_double_layer_potential.hpp"
%include "assembly/helmholtz_3d_hypersingular_operator.hpp"

namespace Bempp
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(Helmholtz3dSingleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(Helmholtz3dDoubleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotential);
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);

} // namespace Bempp

%pythoncode %{

def _constructHelmholtzOperator(className, testSpace, trialSpace, waveNumber):
    basisFunctionType = testSpace.basisFunctionType()
    if (basisFunctionType != trialSpace.basisFunctionType()):
        raise TypeError("BasisFunctionType of testSpace must match that of trialSpace")
    resultType = promoteTypeToComplex(basisFunctionType)
    result = constructObjectTemplatedOnBasis(
        className, basisFunctionType, testSpace, trialSpace, waveNumber)
    result._testSpace = testSpace
    result._trialSpace = trialSpace
    return result

def helmholtz3dSingleLayerPotential(testSpace, trialSpace, waveNumber):
    """Construct a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dSingleLayerPotential", testSpace, trialSpace, waveNumber)

def helmholtz3dDoubleLayerPotential(testSpace, trialSpace, waveNumber):
    """Construct a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dDoubleLayerPotential", testSpace, trialSpace, waveNumber)

def helmholtz3dAdjointDoubleLayerPotential(testSpace, trialSpace, waveNumber):
    """Construct an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dAdjointDoubleLayerPotential", testSpace, trialSpace, waveNumber)

def helmholtz3dHypersingularOperator(testSpace, trialSpace, waveNumber):
    """Construct a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dHypersingularOperator", testSpace, trialSpace, waveNumber)

%}
