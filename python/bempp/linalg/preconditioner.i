%{
#include "linalg/preconditioner.hpp"
%}

namespace Bempp {
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(Preconditioner)
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(acaDiscreteOperatorToPreconditioner)
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_VALUE(acaBlockDiagonalPreconditioner)


%extend Preconditioner {

    %ignore Preconditioner;
    %ignore get;
    %ignore trilinosPreconditioner;
}

BEMPP_EXTEND_CLASS_TEMPLATED_ON_VALUE(Preconditioner)
}

#define shared_ptr boost::shared_ptr
%include "linalg/preconditioner.hpp"
#undef shared_ptr

namespace Bempp {
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(Preconditioner)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(acaDiscreteOperatorToPreconditioner)
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(acaBlockDiagonalPreconditioner)

}

