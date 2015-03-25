#ifndef eigen_types_hpp
#define eigen_types_hpp

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace Bempp {

template <typename T> using Matrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template <typename T> using Vector = Eigen::Matrix<T,Eigen::Dynamic,1>;
template <typename T> using RowVector = Eigen::Matrix<T,1,Eigen::Dynamic>;

// see http://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
// for the following two functions

template <typename T>
void eigenRemoveRowFromMatrix(Matrix<T>& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

template <typename T>
void eigenRemoveColumnFromMatrix(Matrix<T>& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

template <typename T>
Matrix<T> eigenJoinCols(const Matrix<T>& m1, const Matrix<T> m2){

    if (m1.cols() != m2.cols())
        throw std::runtime_error("EigenJoinCols: Column sizes mismatch.");

    Matrix<T> result(m1.rows()+m2.rows(),m1.cols());
    result << m1,m2;
    return result;
}

template <typename T>
Matrix<T> eigenJoinRows(const Matrix<T>& m1, const Matrix<T> m2){

    if (m1.rows() != m2.rows())
        throw std::runtime_error("EigenJoinRows: Row sizes mismatch.");

    Matrix<T> result(m1.rows(),m1.cols()+m2.cols());
    result << m1,m2;
    return result;
}


class EigenInit {

private:

    EigenInit();
    EigenInit(const EigenInit &other);
    const EigenInit &operator=(const EigenInit &other);

    static EigenInit m_singleton;

};

template <typename T>
Matrix<T> eigenMatPinv(const Matrix<T>& mat)
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

