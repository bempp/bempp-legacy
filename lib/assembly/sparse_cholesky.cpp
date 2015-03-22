#include "sparse_cholesky.hpp"

#include "../common/eigen_support.hpp"
#include "../common/boost_make_shared_fwd.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>

namespace Bempp {

shared_ptr<Epetra_CrsMatrix> sparseCholesky(const Epetra_CrsMatrix &mat) {
  // Note: we assume the matrix mat is symmetric and positive-definite
  size_t size = mat.NumGlobalCols();
  if (mat.NumGlobalRows() != size)
    throw std::invalid_argument("sparseCholesky(): matrix must be square");

  int *rowOffsets = 0;
  int *colIndices = 0;
  double *values = 0;
  mat.ExtractCrsDataPointers(rowOffsets, colIndices, values);

  Epetra_SerialComm comm;
  Epetra_LocalMap rowMap(static_cast<int>(size), 0 /* index_base */, comm);
  Epetra_LocalMap columnMap(static_cast<int>(size), 0 /* index_base */, comm);
  shared_ptr<Epetra_CrsMatrix> result = boost::make_shared<Epetra_CrsMatrix>(
      Copy, rowMap, columnMap, mat.GlobalMaxNumEntries());

  Matrix<double> localMat;
  Matrix<double> localCholesky;
  std::vector<bool> processed(size, false);
  for (size_t r = 0; r < size; ++r) {
    if (processed[r])
      continue;
    int localSize = rowOffsets[r + 1] - rowOffsets[r];
    localMat.resize(localSize, localSize);
    localMat.setZero();
    localCholesky.resize(localSize, localSize);
    for (int s = 0; s < localSize; ++s) {
      int row = colIndices[rowOffsets[r] + s];
      for (int c = 0; c < localSize; ++c) {
        int col = colIndices[rowOffsets[row] + c];
        if (col != colIndices[rowOffsets[r] + c])
          throw std::invalid_argument("sparseCholesky(): matrix is not "
                                      "block-diagonal");
        localMat(s, c) = values[rowOffsets[row] + c];
      }
    }
    (localMat - localMat.adjoint()).norm() <
           1e-12 * localMat.norm();
    localCholesky = localMat.llt().matrixL().adjoint();
    for (int s = 0; s < localSize; ++s) {
      int row = colIndices[rowOffsets[r] + s];
      processed[row] = true;
#ifndef NDEBUG
      int errorCode =
#endif
          result->InsertGlobalValues(row, s + 1 /* number of values */,
                                     localCholesky.colptr(s),
                                     colIndices + rowOffsets[r]);
      assert(errorCode == 0);
    }
  }
  result->FillComplete(columnMap, rowMap);

  return result;
}

} // namespace Bempp
