/*// Invocation of this macro is necessary because Helmholtz operators
// have their ResultType determined "implicitly". A better way
// would be to do these things in LinearOperator.__init__,
// but SWIG puts exception-raising code in abstract class constructors...

%define BEMPP_PYTHON_EXTEND_HELMHOLTZ_OPERATOR(CLASS)

%extend CLASS<float>
{
    %pythonappend CLASS
    %{
        self._resultType = "complex64"
    %}
}

%extend CLASS<double>
{
    %pythonappend CLASS
    %{
        self._resultType = "complex128"
    %}
}

%extend CLASS<std::complex<float> >
{
    %pythonappend CLASS
    %{
        self._resultType = "complex64"
    %}
}

%extend CLASS<std::complex<double> >
{
    %pythonappend CLASS
    %{
        self._resultType = "complex128"
    %}
}

%enddef // BEMPP_PYTHON_EXTEND_HELMHOLTZ_OPERATOR*/

%{
#include "assembly/helmholtz_3d_single_layer_potential.hpp"
#include "assembly/helmholtz_3d_double_layer_potential.hpp"
#include "assembly/helmholtz_3d_adjoint_double_layer_potential.hpp"
#include "assembly/helmholtz_3d_hypersingular_operator.hpp"
#include <complex>
%}

// TODO
// %include "helmholtz_3d_operators_docstrings.i"

namespace Bempp
{
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dSingleLayerPotential);
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dDoubleLayerPotential);
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotential);
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);

} // namespace Bempp

%include "assembly/helmholtz_3d_single_layer_potential.hpp"
%include "assembly/helmholtz_3d_double_layer_potential.hpp"
%include "assembly/helmholtz_3d_adjoint_double_layer_potential.hpp"
%include "assembly/helmholtz_3d_hypersingular_operator.hpp"

namespace Bempp
{

/*%extend Helmholtz3dSingleLayerPotential
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
}*/

BEMPP_PYTHON_DECLARE_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dSingleLayerPotential);
BEMPP_PYTHON_DECLARE_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dDoubleLayerPotential);
BEMPP_PYTHON_DECLARE_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotential);
BEMPP_PYTHON_DECLARE_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);

/*BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dSingleLayerPotential);
BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dDoubleLayerPotential);
BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotential);
BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);

BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dSingleLayerPotential);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dDoubleLayerPotential);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotential);
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);*/

} // namespace Bempp

%pythoncode %{

def _constructHelmholtzOperator(className, testSpace, trialSpace, waveNumber):
    basisFunctionType = testSpace._basisFunctionType
    if (basisFunctionType != trialSpace._basisFunctionType):
        raise TypeError("BasisFunctionType of testSpace must match that of trialSpace")
    resultType = promoteTypeToComplex(basisFunctionType)
    return constructObjectTemplatedOnBasis(
        className, basisFunctionType, testSpace, trialSpace, waveNumber)

def helmholtz3dSingleLayerPotential(testSpace, trialSpace, waveNumber, resultType=None):
    """Construct a single-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
        "Helmholtz3dSingleLayerPotential", testSpace, trialSpace, waveNumber)

def helmholtz3dDoubleLayerPotential(testSpace, trialSpace, waveNumber, resultType=None):
    """Construct a double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
    "Helmholtz3dDoubleLayerPotential", testSpace, trialSpace, waveNumber)

def helmholtz3dAdjointDoubleLayerPotential(testSpace, trialSpace, waveNumber, resultType=None):
    """Construct an adjoint double-layer-potential operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
    "Helmholtz3dAdjointDoubleLayerPotential", testSpace, trialSpace, waveNumber)

def helmholtz3dHypersingularOperator(testSpace, trialSpace, waveNumber, resultType=None):
    """Construct a hypersingular operator for the Helmholtz equation in 3D."""
    return _constructHelmholtzOperator(
    "Helmholtz3dHypersingularOperator", testSpace, trialSpace, waveNumber)

%}
