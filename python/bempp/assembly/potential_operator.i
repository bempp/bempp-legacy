%{
#include "assembly/potential_operator.hpp"
%}

BEMPP_DECLARE_SHARED_PTR_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    Bempp::PotentialOperator);

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(PotentialOperator);

#define shared_ptr boost::shared_ptr
%extend PotentialOperator
{

    %apply const arma::Mat<float>& IN_MAT { const arma::Mat<float>& evaluationPoints };
    %apply const arma::Mat<double>& IN_MAT { const arma::Mat<double>& evaluationPoints };
    %apply const shared_ptr<const arma::Mat<float> >& IN_SP_MAT {
        const shared_ptr<const arma::Mat<float> >& evaluationPoints };
    %apply const shared_ptr<const arma::Mat<double> >& IN_SP_MAT {
        const shared_ptr<const arma::Mat<double> >& evaluationPoints };

    %apply arma::Mat<float>& ARGOUT_MAT { arma::Mat<float>& result_ };
    %apply arma::Mat<double>& ARGOUT_MAT { arma::Mat<double>& result_ };
    %apply arma::Mat<std::complex<float> >& ARGOUT_MAT
        { arma::Mat<std::complex<float> >& result_ };
    %apply arma::Mat<std::complex<double> >& ARGOUT_MAT
        { arma::Mat<std::complex<double> >& result_ };

    %ignore evaluateOnGrid;
    %ignore makeEvaluator;

    void _evaluateAtPoints(
        arma::Mat<ResultType>& result_,
        const GridFunction<BasisFunctionType_, ResultType_>& argument,
        const arma::Mat<CoordinateType>& evaluationPoints,
        const Fiber::QuadratureStrategy<
        BasisFunctionType_, ResultType_, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options)
    {
        result_ = $self->evaluateAtPoints(argument, evaluationPoints, quadStrategy, options);
    }

    %ignore evaluateAtPoints;

    %pythoncode {
        def evaluateAtPoints(self, argument, evaluationPoints,
                             evaluationOptions=EvaluationOptions()):
            return self._evaluateAtPoints(argument, evaluationPoints,
                                          self._context.quadStrategy(), evaluationOptions)
    }
}
#undef

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(PotentialOperator);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/potential_operator.hpp"
#undef

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(PotentialOperator);
}

%clear const arma::Mat<float>& evaluationPoints;
%clear const arma::Mat<double>& evaluationPoints;
%clear arma::Mat<float>& result_;
%clear arma::Mat<double>& result_;
%clear arma::Mat<std::complex<float> >& result_;
%clear arma::Mat<std::complex<double> >& result_;
