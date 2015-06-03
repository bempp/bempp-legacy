#ifndef math_helper_impl_hpp
#define math_helper_impl_hpp

#include "math_helper.hpp"

namespace hmat {

template <typename ValueType>
void computeLowRankApproximation(const Matrix<ValueType> &mat, double threshold,
                                 int maxRank, bool &success,
                                 Matrix<ValueType> &U, Matrix<ValueType> &S,
                                 Matrix<ValueType> &V) {
  Eigen::JacobiSVD<Matrix<ValueType>> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // Get the index that is lower than the threshold
  svd.setThreshold(threshold);
  auto rank = svd.rank();

  if (rank > maxRank) {
    success = false;
    return;
  }

  if (rank < std::min(mat.rows(), mat.cols())) {

    U = svd.matrixU().leftCols(rank);
    V = svd.matrixV().leftCols(rank);
    S = svd.singularValues().template cast<ValueType>().head(rank).asDiagonal();

    success = true;

  } else
    success = false;
}

template <typename ValueType>
void compressQB(Matrix<ValueType> &Q, Matrix<ValueType> &B, double threshold, int maxRank,
                bool &success) {

  Matrix<ValueType> U, S, V;

  computeLowRankApproximation(B, threshold, maxRank, success, U, S, V);
  if (success) {
    Q = Q * U * S;
    B = V.adjoint();
  }
}

template <typename ValueType>
void randomizedLowRankApproximation(const matApply_t<ValueType> &applyFun, double threshold,
                                    int maxRank, int sampleDimension, int rows,
                                    bool &success, Matrix<ValueType> &A,
                                    Matrix<ValueType> &B) {

  Matrix<ValueType> Z = Matrix<ValueType>::Random(rows, sampleDimension);
  Eigen::HouseholderQR<Matrix<ValueType>> qr(applyFun(Eigen::Ref<Matrix<ValueType>>(Z), NOTRANS));
  A = qr.householderQ();
  B = applyFun(Eigen::Ref<Matrix<ValueType>>(A), CONJTRANS).adjoint();
  compressQB(A, B, threshold, maxRank, success);
}
}

#endif
