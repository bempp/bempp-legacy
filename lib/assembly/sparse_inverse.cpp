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

namespace Bempp {

shared_ptr<RealSparseMatrix> sparseInverse(const RealSparseMatrix &mat) {
  // Note: we assume the matrix mat is symmetric and positive-definite
  size_t size = mat.cols();
  if (mat.rows() != size)
    throw std::invalid_argument("sparseInverse(): matrix must be square");

  if (!mat.isCompressed())
    throw std::invalid_argument(
        "sparseInverse(): matrix must be in compressed form");

  shared_ptr<RealSparseMatrix> result =
      boost::make_shared<RealSparseMatrix>(size, size);

  std::vector<Eigen::Triplet<double>> triplets;

  auto outerIndexPtr = mat.outerIndexPtr();
  auto innerIndexPtr = mat.innerIndexPtr();

  Matrix<double> localInverse;
  std::vector<bool> processed(size, false);
  for (size_t r = 0; r < size; ++r) {
    if (processed[r])
      continue;
    int localSize = outerIndexPtr[r + 1] - outerIndexPtr[r];
    for (int s = 0; s < localSize; ++s) {
      int col = innerIndexPtr[outerIndexPtr[r] + s];
      for (int c = 0; c < localSize; ++c) {
        int row = innerIndexPtr[outerIndexPtr[col] + c];
        if (row != innerIndexPtr[outerIndexPtr[r] + c])
          throw std::invalid_argument(
              "sparseInverse(): matrix is not block-diagonal. "
              "If this error occurs during the assembly of a "
              "synthetic boundary operator, make sure that you "
              "set the internalTrialSpace and internalTestSpace "
              "parameters in its constructor to a "
              "\"discontinuous\" function space, i.e. with each "
              "of its basis functions living on a single "
              "element only");
      }
    }
    localInverse =
        Matrix<double>(mat.block(r, r, localSize, localSize)).inverse();
    for (int s = 0; s < localSize; ++s) {
      processed[r + s] = true;
      for (int c = 0; c < localSize; ++c)
        triplets.push_back(
            Eigen::Triplet<double>(r + s, r + c, localInverse(s, c)));
    }
  }
  result->setFromTriplets(triplets.begin(), triplets.end());
  return result;
}

} // namespace Bempp
