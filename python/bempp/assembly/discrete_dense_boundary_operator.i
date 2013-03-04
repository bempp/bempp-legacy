%{
#include "assembly/discrete_dense_boundary_operator.hpp"
#include "fiber/scalar_traits.hpp"
#include "common/shared_ptr.hpp"
%}

namespace Bempp 
{
    %ignore DiscreteDenseBoundaryOperator;

    %apply const arma::Mat<float>& IN_MAT {
        const arma::Mat<float>& mat
    };
    %apply const arma::Mat<double>& IN_MAT {
        const arma::Mat<double>& mat
    };
    %apply const arma::Mat<std::complex<float> >& IN_MAT {
        const arma::Mat<std::complex<float> >& mat
    };
    %apply const arma::Mat<std::complex<double> >& IN_MAT {
        const arma::Mat<std::complex<double> >& mat
    };
}

#define shared_ptr boost::shared_ptr
//%include "assembly/discrete_dense_boundary_operator.hpp"
namespace Bempp 
{
template <typename ValueType>
shared_ptr<DiscreteBoundaryOperator<ValueType> > discreteDenseBoundaryOperator(
        const arma::Mat<ValueType>& mat);
}
#undef boost::shared_ptr

namespace Bempp 
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(
    discreteDenseBoundaryOperator);

%clear const arma::Mat<float>& mat;
%clear const arma::Mat<double>& mat;
%clear const arma::Mat<std::complex<float> >& mat;
%clear const arma::Mat<std::complex<float> >& mat;

}
