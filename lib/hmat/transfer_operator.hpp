// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_TRANSFER_OPERATOR
#define HMAT_TRANSFER_OPERATOR

#include "common.hpp"
#include "scalar_traits.hpp"
#include "eigen_fwd.hpp"

#include <vector>

namespace hmat {

template <typename ValueType> class TransferOperator {
public:
  virtual void apply(const Eigen::Ref<Matrix<ValueType>> &X,
                     Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
                     ValueType alpha, ValueType beta) const = 0;

};
}

#endif
