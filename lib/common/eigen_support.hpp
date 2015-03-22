#ifndef eigen_types_hpp
#define eigen_types_hpp

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace Bempp {

template <typename T> using Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template <typename T> using Vector = Eigen::Matrix<T,Eigen::Dynamic,1>;

class EigenInit {

private:

    EigenInit();
    EigenInit(const EigenInit &other);
    const EigenInit &operator=(const EigenInit &other);

    static EigenInit m_singleton;

};

template <typename T>
Matrix<T> eigenMatPinv(const Matrix<T>& mat) const
{

    if (mat.rows()==mat.cols())
        return mat.inverse();

    if (mat.rows()>mat.cols())
        return (mat.transpose()*mat).inverse()*mat.transpose();
    else
        return mat.transpose()*(mat*mat.transpose()).inverse();
}


}


#endif // eigen_types_hpp

