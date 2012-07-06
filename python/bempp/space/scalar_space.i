%{
#include "space/space.hpp"
#include <complex>
%}

%shared_ptr(Bempp::ScalarSpace<float>);
%shared_ptr(Bempp::ScalarSpace<double>);
%shared_ptr(Bempp::ScalarSpace<std::complex<float> >);
%shared_ptr(Bempp::ScalarSpace<std::complex<double> >);

// TODO
// %include "grid_docstrings.i"


namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(ScalarSpace);
}

%include "space/scalar_space.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(ScalarSpace);
}
