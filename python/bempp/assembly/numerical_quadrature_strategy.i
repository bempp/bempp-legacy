%{
#include "assembly/numerical_quadrature_strategy.hpp"
%}

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    NumericalQuadratureStrategy);
BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    NumericalQuadratureStrategy);
}

%include "assembly/numerical_quadrature_strategy.hpp"

%shared_ptr(Bempp::NumericalQuadratureStrategy<float, float>);
%shared_ptr(Bempp::NumericalQuadratureStrategy<float, std::complex<float> >);
%shared_ptr(Bempp::NumericalQuadratureStrategy<std::complex<float>, std::complex<float> >);
%shared_ptr(Bempp::NumericalQuadratureStrategy<double, double>);
%shared_ptr(Bempp::NumericalQuadratureStrategy<double, std::complex<double> >);
%shared_ptr(Bempp::NumericalQuadratureStrategy<std::complex<double>, std::complex<double> >);

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(
NumericalQuadratureStrategy);
}

%pythoncode %{

def numericalQuadratureStrategy(basisFunctionType, resultType, accuracyOptions=None):
    """Construct a NumericalQuadratureStrategy object.

    *Arguments:*
        basisFunctionType (string)
            Type used to represent values of basis functions.

        resultType (string)
            Type used to represent values of boundary-element integrals.

        accuracyOptions (AccuracyOptions)
            Determines quadrature order. If set to None, default quadrature orders
            are used.

        The following combinations of basisFunctionType and resultType are allowed:

            basisFunctionType     resultType
            --------------------------------
            "float32"             "float32" or "complex64"
            "float64"             "float64" or "complex128"
            "complex64"           "complex64"
            "complex128"          "complex128"
    """
    basisFunctionType = checkType(basisFunctionType)
    resultType = checkType(resultType)
    if accuracyOptions is None:
        accuracyOptions = AccuracyOptions()
    name = 'NumericalQuadratureStrategy'
    return constructObjectTemplatedOnBasisAndResult(
        name, basisFunctionType, resultType, accuracyOptions)

%}
