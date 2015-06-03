#ifndef math_helper_hpp
#define math_helper_hpp

#include "common.hpp"
#include "eigen_fwd.hpp"
#include <functional>

namespace hmat {

template <typename ValueType>
using matApply_t = std::function<
    Matrix<ValueType>(const Eigen::Ref<Matrix<ValueType>> &, const TransposeMode trans)>;

template <typename ValueType>
void computeLowRankApproximation(const Matrix<ValueType> &mat, double threshold,
                                 int maxRank, bool &accepted,
                                 Matrix<ValueType> &U, Matrix<ValueType> &S,
                                 Matrix<ValueType> &V);

template <typename ValueType>
void compressQB(Matrix<ValueType> &Q, Matrix<ValueType> &B, double threshold,
                int maxRank, bool &success);

template <typename ValueType>
void randomizedLowRankApproximation(const matApply_t<ValueType> &applyFun, double threshold,
                                    int maxRank, int sampleDimension, int rows,
                                    bool &success, Matrix<ValueType> &A,
                                    Matrix<ValueType> &B);
}

#include "math_helper_impl.hpp"

#endif
