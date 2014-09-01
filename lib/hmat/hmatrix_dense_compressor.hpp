// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DENSE_COMPRESSOR_HPP
#define HMAT_HMATRIX_DENSE_COMPRESSOR_HPP

#include "common.hpp"
#include "hmatrix_compressor.hpp"
#include "data_accessor.hpp"

namespace hmat {

template <typename ValueType, int N>
class HMatrixDenseCompressor : public HMatrixCompressor<ValueType, N> {
public:
  HMatrixDenseCompressor(const DataAccessor<ValueType, N> &dataAccessor);

  void compressBlock(const BlockClusterTreeNode<N> &blockClusterTreeNode,
                     shared_ptr<HMatrixData<ValueType>> &hMatrixData) const
      override;

private:
  const DataAccessor<ValueType, N> &m_dataAccessor;
};
}

#include "hmatrix_dense_compressor_impl.hpp"

#endif
