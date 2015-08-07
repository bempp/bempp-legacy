// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_LEAF_BASIS
#define HMAT_LEAF_BASIS

#include "common.hpp"
#include "scalar_traits.hpp"
#include "eigen_fwd.hpp"

#include <vector>

namespace hmat {

template <typename ValueType> class LeafBasis {
public:
  virtual void apply(const Eigen::Ref<Matrix<ValueType>> &X,
                     Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
                     ValueType alpha, ValueType beta) const = 0;

  virtual int rank() const;

};
}

#endif
