%{
#include "assembly/elementary_potential_operator.hpp"
%}


namespace Bempp 
{

template<typename BasisFunctionType_, typename KernelType_, typename ResultType_> class ElementaryPotentialOperator;

BEMPP_EXTEND_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryPotentialOperator);

}

%include "assembly/elementary_potential_operator.hpp"

namespace Bempp 
{

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryPotentialOperator);

}
