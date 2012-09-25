%{
#include "assembly/discrete_inverse_sparse_boundary_operator.hpp"
  %}

namespace Bempp {

#define shared_ptr boost::shared_ptr
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> >
discreteSparseInverse(const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& discreteOp);
#undef shared_ptr

BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_VALUE(discreteSparseInverse)

}
