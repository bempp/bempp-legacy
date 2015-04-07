#include <Eigen/Dense>

namespace hmat {

template <typename T> using Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template <typename T> using Vector = Eigen::Matrix<T,Eigen::Dynamic,1>;
template <typename T> using RowVector = Eigen::Matrix<T,1,Eigen::Dynamic>;

}
