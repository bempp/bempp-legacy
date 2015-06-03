#ifndef math_helper_impl_hpp
#define math_helper_impl_hpp

#include "math_helper.hpp"

namespace hmat {


    template <typename ValueType>
    void computeLowRankApproximation(
        const Matrix<ValueType>& mat, double threshold, 
        int maxRank, bool& success,
        Matrix<ValueType>& U, Matrix<ValueType>& S, Matrix<ValueType>& V)
    {
      Eigen::JacobiSVD svd(
              mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
      // Get the index that is lower than the threshold
      svd.setThreshold(threshold);
      auto rank = rank();

      if (rank > maxRank) {
          success = false;
          return;
      }

      if (rank<std::min(mat.rows(),mat.cols())){

          U = svd.matrixU().leftCols(rank);
          V = svd.matrixV().leftCols(rank);
          S = svd.singularValues().head(rank).asDiagonal();

          success = true;

      }
      else
          success = false;

    }

    template <typename ValueType>
    void compressQB(Matrix<ValueType>& Q, Matrix<ValueType>& B,
            double threshold, bool& success){


        Matrix<ValueType> U, S, V;

        computeLowRankApproximation(B,threshold,success,U,S,V);
        if (success){
            Q = Q*U*S;
            B = V.adjoint();
        }

    }
            
    template <typename ValueType>
    void randomizedLowRankApproximation(const matApply& applyFun,
            double threshold, int maxRank, int sampleDimension, 
            int rows, bool& success, 
            Matrix<ValueType>& A, Matrix<ValueType>& B){

        Matrix<ValueType> Z = Matrix<ValueType>::Random(rows,
                sampleDimension);
        Eigen::HouseholderQR qr(applyFun(Z,NOTRANS));
        Matrix<ValueType> A = qr.householderQ();
        B = applyFun(Q,CONJTRANS).adjoint();
        compressQB(A,B,threshold,success);

    }

}




#endif
