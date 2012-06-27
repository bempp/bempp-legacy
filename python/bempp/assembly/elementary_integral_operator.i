%{
#include "assembly/elementary_integral_operator.hpp"
%}

// TODO
// %include "elementary_integral_operator_docstrings.i"

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
ElementaryIntegralOperator);

%extend ElementaryIntegralOperator {

    %apply const arma::Mat<float>& IN_MAT {
        const arma::Mat<float>& evaluationPoints
    };
    %apply const arma::Mat<double>& IN_MAT {
        const arma::Mat<double>& evaluationPoints
    };

    %apply arma::Mat<float>& ARGOUT_MAT {
        arma::Mat<float>& result_
    };
    %apply arma::Mat<double>& ARGOUT_MAT {
        arma::Mat<double>& result_
    };
    %apply arma::Mat<std::complex<float> >& ARGOUT_MAT {
        arma::Mat<std::complex<float> >& result_
    };
    %apply arma::Mat<std::complex<double> >& ARGOUT_MAT {
        arma::Mat<std::complex<double> >& result_
    };

    void applyOffSurface(
        arma::Mat<ResultType>& result_,
        const GridFunction<BasisFunctionType, ResultType>& argument,
        const arma::Mat<CoordinateType>& evaluationPoints,
        const LocalAssemblerFactory& factory,
        const EvaluationOptions& options = EvaluationOptions()) const
    {
        result_ = $self->applyOffSurface(argument, evaluationPoints, factory, options);
    }

    %ignore applyOffSurface;
}

} // namespace Bempp

%include "assembly/elementary_integral_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
    ElementaryIntegralOperator);

%clear arma::Mat<float>& result_;
%clear arma::Mat<double>& result_;
%clear arma::Mat<std::complex<float> >& result_;
%clear arma::Mat<std::complex<double> >& result_;

%clear const arma::Mat<float>& evaluationPoints;
%clear const arma::Mat<double>& evaluationPoints;

}
