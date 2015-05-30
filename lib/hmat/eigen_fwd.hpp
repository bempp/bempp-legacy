#ifndef eigen_fwd_hpp
#define eigen_fwd_hpp

#include <Eigen/Dense>

namespace hmat {

template <typename T> using Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template <typename T> using Vector = Eigen::Matrix<T,Eigen::Dynamic,1>;
template <typename T> using RowVector = Eigen::Matrix<T,1,Eigen::Dynamic>;

template <typename Derived>
bool isnan(const Eigen::MatrixBase<Derived>& matrix)
{
    return (!((matrix-matrix).array()==(matrix-matrix).array()).all());

}

}

#endif
