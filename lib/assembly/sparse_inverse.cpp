// Copyright (C) 2011-2013 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "sparse_inverse.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/eigen_support.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_LocalMap.h>
#include <Epetra_SerialComm.h>

namespace Bempp {

shared_ptr<Epetra_CrsMatrix> sparseInverse(const Epetra_CrsMatrix &mat) {
  // Note: we assume the matrix mat is symmetric and positive-definite
  size_t size = mat.NumGlobalCols();
  if (mat.NumGlobalRows() != size)
    throw std::invalid_argument("sparseInverse(): matrix must be square");

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
  Matrix<double> localInverse;
  std::vector<bool> processed(size, false);
  for (size_t r = 0; r < size; ++r) {
    if (processed[r])
      continue;
    int localSize = rowOffsets[r + 1] - rowOffsets[r];
    localMat.resize(localSize, localSize);
    localMat.setZero();
    localInverse.resize(localSize, localSize);
    for (int s = 0; s < localSize; ++s) {
      int row = colIndices[rowOffsets[r] + s];
      for (int c = 0; c < localSize; ++c) {
        int col = colIndices[rowOffsets[row] + c];
        if (col != colIndices[rowOffsets[r] + c])
          throw std::invalid_argument(
              "sparseInverse(): matrix is not block-diagonal. "
              "If this error occurs during the assembly of a "
              "synthetic boundary operator, make sure that you "
              "set the internalTrialSpace and internalTestSpace "
              "parameters in its constructor to a "
              "\"discontinuous\" function space, i.e. with each "
              "of its basis functions living on a single "
              "element only");
        localMat(s, c) = values[rowOffsets[row] + c];
      }
    }
    localInverse = localMat.inverse();
    for (int s = 0; s < localSize; ++s) {
      int row = colIndices[rowOffsets[r] + s];
      processed[row] = true;
#ifndef NDEBUG
      int errorCode =
#endif
          result->InsertGlobalValues(row, localSize /* number of values */,
                                     localInverse.colptr(s),
                                     colIndices + rowOffsets[r]);
      assert(errorCode == 0);
    }
  }
  result->FillComplete(columnMap, rowMap);

  return result;
}

} // namespace Bempp
