#ifndef HMAT_MATH_HELPER_HPP
#define HMAT_MATH_HELPER_HPP

#include "common.hpp"
#include "eigen_fwd.hpp"
#include <functional>

namespace hmat {

template <typename ValueType>
using matApply_t = std::function<Matrix<ValueType>(
    const Eigen::Ref<Matrix<ValueType>> &, const TransposeMode trans)>;

template <typename ValueType>
void computeLowRankApproximation(const Matrix<ValueType> &mat, double threshold,
                                 int maxRank, bool &accepted,
                                 Matrix<ValueType> &U, Matrix<ValueType> &S,
                                 Matrix<ValueType> &V);

template <typename ValueType>
void compressQB(Matrix<ValueType> &Q, Matrix<ValueType> &B, double threshold,
                int maxRank, bool &success);

template <typename ValueType>
void randomizedLowRankApproximation(const matApply_t<ValueType> &applyFun,
                                    int rows, int cols, double threshold,
                                    int maxRank, int sampleDimension,
                                    bool &success, Matrix<ValueType> &A,
                                    Matrix<ValueType> &B);

template <typename ValueType>
std::size_t computeRank(const Eigen::JacobiSVD<Matrix<ValueType>> &svd,
                        double threshold);
}

#include "math_helper_impl.hpp"

#endif
