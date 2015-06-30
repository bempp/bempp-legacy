// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_ACA_COMPRESSOR_HPP
#define HMAT_HMATRIX_ACA_COMPRESSOR_HPP

#include "common.hpp"
#include "eigen_fwd.hpp"
#include "hmatrix_compressor.hpp"
#include "hmatrix_dense_compressor.hpp"
#include "data_accessor.hpp"
#include <set>

namespace hmat {

template <typename ValueType, int N>
class HMatrixAcaCompressor : public HMatrixCompressor<ValueType, N> {
public:
  HMatrixAcaCompressor(const DataAccessor<ValueType, N> &dataAccessor,
                       double eps, unsigned int maxRank);

  void
  compressBlock(const BlockClusterTreeNode<N> &blockClusterTreeNode,
                shared_ptr<HMatrixData<ValueType>> &hMatrixData) const override;

  enum class ModeType { ROW, COL };
  enum class CrossStatusType { SUCCESS, ZERO };
  enum class AcaStatusType { SUCCESS, ZERO_TERMINATION, RANK_LIMIT_REACHED };

private:
  static std::size_t randomIndex(const IndexRangeType &range,
                                 std::set<std::size_t> &previousIndices);

  CrossStatusType computeCross(
      const BlockClusterTreeNode<N> &blockClusterTreeNode,
      const Matrix<ValueType> &A, const Matrix<ValueType> &B,
      std::size_t &nextPivot, Matrix<ValueType> &origRow,
      Matrix<ValueType> &origCol, Matrix<ValueType> &row,
      Matrix<ValueType> &col, std::vector<std::size_t> &rowApproxCounter,
      std::vector<std::size_t> &colApproxCounter, std::vector<double> &rowNorms,
      std::vector<double> &colNorms, ModeType mode, double zeroTol) const;

  bool selectPivotWithMinNorm2(const std::vector<double> &norm2,
                               const std::vector<int> &apprx_times, int &pivot);

  double updateLowRankBlocksAndNorm(const Matrix<ValueType> &newRow,
                                    const Matrix<ValueType> &newCol,
                                    Matrix<ValueType> &A, Matrix<ValueType> &B,
                                    double currentBlockNorm) const;

  bool checkConvergence(const Matrix<ValueType> &row,
                        const Matrix<ValueType> &col, double tol) const;

  AcaStatusType aca(const BlockClusterTreeNode<N> &blockClusterTreeNode,
                    std::size_t startPivot, Matrix<ValueType> &A,
                    Matrix<ValueType> &B, std::size_t maxIterations,
                    std::vector<size_t> &rowApproxCounter,
                    std::vector<size_t> &colApproxCounter,
                    std::vector<double> &rowNorms,
                    std::vector<double> &colNorms, double &blockNorm,
                    double eps, double zeroTol, ModeType mode) const;

  const DataAccessor<ValueType, N> &m_dataAccessor;
  double m_eps;
  unsigned int m_maxRank;
  HMatrixDenseCompressor<ValueType, N> m_hMatrixDenseCompressor;
};
}

#include "hmatrix_aca_compressor_impl.hpp"

#endif
