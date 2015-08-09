#ifndef HMAT_MATH_HELPER_IMPL_HPP
#define HMAT_MATH_HELPER_IMPL_HPP

#include "math_helper.hpp"
#include <limits>

namespace hmat {

template <typename ValueType>
void computeLowRankApproximation(const Matrix<ValueType> &mat, double threshold,
                                 int maxRank, bool &success,
                                 Matrix<ValueType> &U, Matrix<ValueType> &S,
                                 Matrix<ValueType> &V) {
  if (!mat.allFinite())
    std::cout << "mat is infinite" << std::endl;
  Eigen::JacobiSVD<Matrix<ValueType>> svd(mat, Eigen::ComputeThinU |
                                                   Eigen::ComputeThinV);
  // Get the index that is lower than the threshold
  auto rank = computeRank(svd,threshold);
  if (rank == 0)
    rank = 1;

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
void compressQB(Matrix<ValueType> &Q, Matrix<ValueType> &B, double threshold,
                int maxRank, bool &success) {

  Matrix<ValueType> U, S, V;

  computeLowRankApproximation(B, threshold, maxRank, success, U, S, V);
  if (success) {
    Q = Q * U * S;
    B = V.adjoint();
  }
}

template <typename ValueType>
void randomizedLowRankApproximation(const matApply_t<ValueType> &applyFun,
                                    int rows, int cols, double threshold,
                                    int maxRank, int sampleDimension,
                                    bool &success, Matrix<ValueType> &A,
                                    Matrix<ValueType> &B) {

  int actual_sample_size = std::min(rows, sampleDimension);
  Matrix<ValueType> identity =
      Matrix<ValueType>::Identity(rows, actual_sample_size);
  Matrix<ValueType> Z = Matrix<ValueType>::Random(cols, actual_sample_size);

  A = applyFun(Eigen::Ref<Matrix<ValueType>>(Z),NOTRANS).
      colPivHouseholderQr().householderQ()*identity;

  B = applyFun(Eigen::Ref<Matrix<ValueType>>(A), CONJTRANS).adjoint();
  compressQB(A, B, threshold, maxRank, success);
}


template <typename ValueType>
std::size_t computeRank(const Eigen::JacobiSVD<Matrix<ValueType>>& svd, double threshold){

    auto singularValues = svd.singularValues().array();
    auto maxSingVal = singularValues(0);
    return (singularValues>(threshold*maxSingVal)).count();


}

}

#endif
