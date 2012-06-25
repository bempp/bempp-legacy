%{
#include "fiber/scalar_traits.hpp"
%}

%include "fiber/scalar_traits.hpp"
%include "common/scalar_traits.hpp"

namespace Fiber
{

%template() ScalarTraits<float>;
%template() ScalarTraits<double>;
%template() ScalarTraits<std::complex<float> >;
%template() ScalarTraits<std::complex<double> >;

} // namespace Fiber
