%{
#include "assembly/discrete_boundary_operator.hpp"
#include <complex>
%}

// TODO
// %include "discrete_boundary_operator_docstrings.i"

%shared_ptr(Thyra::LinearOpDefaultBase<float>);
%shared_ptr(Thyra::LinearOpDefaultBase<double>);
%shared_ptr(Thyra::LinearOpDefaultBase<std::complex<float> >);
%shared_ptr(Thyra::LinearOpDefaultBase<std::complex<double> >);

%shared_ptr(Bempp::DiscreteBoundaryOperator<float>);
%shared_ptr(Bempp::DiscreteBoundaryOperator<double>);
%shared_ptr(Bempp::DiscreteBoundaryOperator<std::complex<float> >);
%shared_ptr(Bempp::DiscreteBoundaryOperator<std::complex<double> >);

namespace Thyra
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(LinearOpDefaultBase);

template <typename ValueType>
class LinearOpDefaultBase
{
public:
    virtual ~LinearOpDefaultBase() = 0; // prevent instantiation
};

} // namespace Thyra


namespace Bempp
{

DECLARE_TEMPLATE_VALUE_METHOD_AUTO_DOCSTRING(DiscreteBoundaryOperator, apply);

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator);

%extend DiscreteBoundaryOperator
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
    %apply const arma::Col<std::complex<double> >& IN_COL {
        const arma::Col<std::complex<double> >& x_in
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

%include "assembly/discrete_boundary_operator.hpp"

namespace Thyra
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(LinearOpDefaultBase);
}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(DiscreteBoundaryOperator);

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
%clear arma::Mat<std::complex<double> >& mat_out;

} // namespace Bempp

