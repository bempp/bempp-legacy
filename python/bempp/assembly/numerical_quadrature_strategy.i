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

def numericalQuadratureStrategy(accuracyOptions=None,
        basisFunctionType='float64',resultType='float64'):
    if basisFunctionType is not None:
        basisFunctionType = checkType(basisFunctionType)
    if resultType is not None:
        resultType = checkType(resultType)
    if accuracyOptions is None:
        accuracyOptions = AccuracyOptions()
    name = 'NumericalQuadratureStrategy'
    return constructObjectTemplatedOnBasisAndResult(
        name, basisFunctionType, resultType, accuracyOptions)

%}
