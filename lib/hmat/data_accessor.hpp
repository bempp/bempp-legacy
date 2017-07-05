// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_DATA_ACCESSOR_HPP
#define HMAT_DATA_ACCESSOR_HPP

#include "common.hpp"
#include "eigen_fwd.hpp"

namespace hmat {

template <typename ValueType, int N>
class DataAccessor {
  public:
  virtual void
  computeMatrixBlock(const IndexRangeType& rowIndexRange,
      const IndexRangeType& columnIndexRange,
      const BlockClusterTreeNode<N>& blockClusterTreeNode,
      Matrix<ValueType>& data) const = 0;

  virtual double scale(const BlockClusterTreeNode<N>& node) const = 0;

  virtual void dofVolumes(Vector<double>& testVolumes, Vector<double>& trialVolumes) const = 0;
};
}

#endif
