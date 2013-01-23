%{
#include "assembly/assembled_potential_operator.hpp"
%}

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AssembledPotentialOperator);

%extend AssembledPotentialOperator
{
    %apply arma::Mat<float>& ARGOUT_MAT { arma::Mat<float>& result_ };
    %apply arma::Mat<double>& ARGOUT_MAT { arma::Mat<double>& result_ };
    %apply arma::Mat<std::complex<float> >& ARGOUT_MAT
        { arma::Mat<std::complex<float> >& result_ };
    %apply arma::Mat<std::complex<double> >& ARGOUT_MAT
        { arma::Mat<std::complex<double> >& result_ };

    void _apply(
        arma::Mat<ResultType>& result_,
        const GridFunction<BasisFunctionType, ResultType>& argument)
    {
        result_ = $self->apply(argument);
    }

    %ignore apply;

    %pythoncode {
        def apply(self, argument):
            return self._apply(argument)
    }
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AssembledPotentialOperator);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/assembled_potential_operator.hpp"
#undef

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(AssembledPotentialOperator);
}

%clear arma::Mat<float>& result_;
%clear arma::Mat<double>& result_;
%clear arma::Mat<std::complex<float> >& result_;
%clear arma::Mat<std::complex<double> >& result_;
