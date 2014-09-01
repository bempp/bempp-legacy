// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_ACA_COMPRESSOR_HPP
#define HMAT_HMATRIX_ACA_COMPRESSOR_HPP

#include "common.hpp"
#include "hmatrix_compressor.hpp"
#include "hmatrix_dense_compressor.hpp"
#include "data_accessor.hpp"

namespace hmat {

template <typename ValueType, int N>
class HMatrixAcaCompressor : public HMatrixCompressor<ValueType, N> {
public:
  HMatrixAcaCompressor(const DataAccessor<ValueType, N> &dataAccessor,
                       double eps, unsigned int maxRank);

  void compressBlock(const BlockClusterTreeNode<N> &blockClusterTreeNode,
                     shared_ptr<HMatrixData<ValueType>> &hMatrixData) const
      override;

private:
  void evaluateMatMinusLowRank(const BlockClusterTreeNode<N> &blockClusterTreeNode,
      const IndexRangeType& rowIndexRange,
      const IndexRangeType& columnIndexRange,
      arma::Mat<ValueType> &data, 
      const std::vector<shared_ptr<arma::Mat<ValueType>>> &previousColumns,
      const std::vector<shared_ptr<arma::Mat<ValueType>>> &previousRows) const;

  static std::size_t intRand(const IndexRangeType& range);

  const DataAccessor<ValueType, N> &m_dataAccessor;
  double m_eps;
  unsigned int m_maxRank;
  HMatrixDenseCompressor<ValueType,N> m_hMatrixDenseCompressor;
};
}

#include "hmatrix_aca_compressor_impl.hpp"

#endif
