%{
#include "assembly/scaled_abstract_boundary_operator.hpp"
%}

// TODO
// %include "scaled_abstract_boundary_operator_docstrings.i"

%include "assembly/scaled_abstract_boundary_operator.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_AND_RESULT(ScaledAbstractBoundaryOperator);
}
