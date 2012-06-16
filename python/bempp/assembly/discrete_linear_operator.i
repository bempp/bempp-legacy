%{
#include "assembly/discrete_linear_operator.hpp"
#include <complex>
%}

// TODO
// %include "discrete_linear_operator_docstrings.i"

namespace Bempp
{

DECLARE_TEMPLATE_VALUE_METHOD_AUTO_DOCSTRING(DiscreteLinearOperator, apply);

BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(DiscreteLinearOperator);

%extend DiscreteLinearOperator
{
    // this function is only for internal use
    %ignore addBlock;

    %apply const arma::Col<float>& IN_COL {
        const arma::Col<float>& x_in
    };
    %apply const arma::Col<double>& IN_COL {
        const arma::Col<double>& x_in
    };
    %apply const arma::Col<std::complex<float> >& IN_COL {
        const arma::Col<std::complex<float> >& x_in
    };
    %apply const arma::Col<std::complex<double>& IN_COL {
        const arma::Col<std::complex<double>& x_in
    };

    %apply arma::Col<float>& INPLACE_COL {
        arma::Col<float>& y_inout
    };
    %apply arma::Col<double>& INPLACE_COL {
        arma::Col<double>& y_inout
    };
    %apply arma::Col<std::complex<float> >& INPLACE_COL {
        arma::Col<std::complex<float> >& y_inout
    };
    %apply arma::Col<std::complex<double> >& INPLACE_COL {
        arma::Col<std::complex<double> >& y_inout
    };

    %apply arma::Mat<float>& ARGOUT_MAT {
        arma::Mat<float>& mat_out
    };
    %apply arma::Mat<double>& ARGOUT_MAT {
        arma::Mat<double>& mat_out
    };
    %apply arma::Mat<std::complex<float> >& ARGOUT_MAT {
        arma::Mat<std::complex<float> >& mat_out
    };
    %apply arma::Mat<std::complex<double> >& ARGOUT_MAT {
        arma::Mat<std::complex<double> >& mat_out
    };
    void asMatrix(arma::Mat<ValueType>& mat_out)
    {
        mat_out = $self->asMatrix();
    }

    %ignore asMatrix;
}

} // namespace Bempp

%include "assembly/discrete_linear_operator.hpp"

namespace Bempp
{
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_VALUE(DiscreteLinearOperator);

%clear const arma::Col<float>& x_in;
%clear const arma::Col<double>& x_in;
%clear const arma::Col<std::complex<float> >& x_in;
%clear const arma::Col<std::complex<double> >& x_in;

%clear arma::Col<float>& y_inout;
%clear arma::Col<double>& y_inout;
%clear arma::Col<std::complex<float> >& y_inout;
%clear arma::Col<std::complex<double> >& y_inout;

%clear arma::Mat<float>& mat_out;
%clear arma::Mat<double>& mat_out;
%clear arma::Mat<std::complex<float> >& mat_out;
%clear arma::Mat<std::complex<float> >& mat_out;

/*%clear const arma::Col<ValueType>& x_in;
%clear arma::Col<ValueType>& y_inout;*/
} // namespace Bempp

