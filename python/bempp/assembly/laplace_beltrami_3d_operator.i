%{
#include "assembly/laplace_beltrami_3d_operator.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") laplaceBeltrami3dOperator;

// Do not emit warnings about ignored assignment operators.
%warnfilter(362) LaplaceBeltrami3dOperator::operator=;
}

#define shared_ptr boost::shared_ptr
%include "assembly/laplace_beltrami_3d_operator.hpp"
#undef shared_ptr

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(laplaceBeltrami3dOperator);
}
