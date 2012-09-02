%{
#include "assembly/context.hpp"
%}

// TODO
// %include "context_docstrings.i"

%shared_ptr(Bempp::Context<float, float>);
%shared_ptr(Bempp::Context<float, std::complex<float> >);
%shared_ptr(Bempp::Context<std::complex<float>, std::complex<float> >);
%shared_ptr(Bempp::Context<double, double>);
%shared_ptr(Bempp::Context<double, std::complex<double> >);
%shared_ptr(Bempp::Context<std::complex<double>, std::complex<double> >);

namespace Bempp
{

BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

class AssemblyOptions;

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Context);

} // namespace Bempp

#define shared_ptr boost::shared_ptr
%include "assembly/context.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(Context);
}
