// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "bempp/common/config_trilinos.hpp"
#include "bempp/common/config_ahmed.hpp"

#include "discrete_sparse_boundary_operator.hpp"

#include "ahmed_mblock_array_deleter.hpp"
#include "discrete_aca_boundary_operator.hpp"
#include "index_permutation.hpp"
#include "sparse_to_h_matrix_converter.hpp"

#include "../common/boost_shared_array_fwd.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/parallelization_options.hpp"

#include <iostream>
#include <stdexcept>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialComm.h>

namespace Bempp {

// Helper functions for the applyBuiltIn member function
namespace {

template <typename ValueType>
void reallyApplyBuiltInImpl(const Epetra_CrsMatrix &mat,
                            const TranspositionMode trans,
                            const Eigen::Ref<Vector<ValueType>> &x_in,
                            Eigen::Ref<Vector<ValueType>> y_inout,
                            const ValueType alpha, const ValueType beta);

template <>
void reallyApplyBuiltInImpl<double>(const Epetra_CrsMatrix &mat,
                                    const TranspositionMode trans,
                                    const Eigen::Ref<Vector<double>> &x_in,
                                    Eigen::Ref<Vector<double>> y_inout,
                                    const double alpha, const double beta) {
  if (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE) {
    assert(mat.NumGlobalRows() == static_cast<int>(x_in.rows()));
    assert(mat.NumGlobalCols() == static_cast<int>(y_inout.rows()));
  } else {
    assert(mat.NumGlobalCols() == static_cast<int>(x_in.rows()));
    assert(mat.NumGlobalRows() == static_cast<int>(y_inout.rows()));
  }

  Epetra_Map map_x((int)x_in.rows(), 0, Epetra_SerialComm());
  Epetra_Map map_y((int)y_inout.rows(), 0, Epetra_SerialComm());

  Epetra_Vector vec_x(View, map_x, const_cast<double *>(x_in.data()));
  // vec_temp will store the result of matrix * x_in
  Epetra_Vector vec_temp(map_y, false /* no need to initialise to zero */);

  mat.Multiply(trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE, vec_x,
               vec_temp);

  if (beta == 0.)
    for (size_t i = 0; i < y_inout.rows(); ++i)
      y_inout(i) = alpha * vec_temp[i];
  else
    for (size_t i = 0; i < y_inout.rows(); ++i)
      y_inout(i) = alpha * vec_temp[i] + beta * y_inout(i);
}

template <>
void reallyApplyBuiltInImpl<float>(const Epetra_CrsMatrix &mat,
                                   const TranspositionMode trans,
                                   const Eigen::Ref<Vector<float>> &x_in,
                                   Eigen::Ref<Vector<float>> y_inout, const float alpha,
                                   const float beta) {
  // Copy the float vectors to double vectors
  Vector<double> x_in_double(x_in.rows());
  for (int i = 0; i < x_in.rows(); ++i) x_in_double(i) = x_in(i);

  Vector<double> y_inout_double(y_inout.rows());
  if (beta != 0.f)
      for (int i = 0 ;i < y_inout.rows(); ++i ) y_inout_double(i) = y_inout(i);

  // Do the operation on the double vectors
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_in_double), Eigen::Ref<Vector<double>>(y_inout_double), alpha,
                                 beta);

  // Copy the result back to the float vector
  for (int i = 0; i < y_inout_double.rows(); ++i) y_inout(i) = y_inout_double(i);

}

template <>
void reallyApplyBuiltInImpl<std::complex<float>>(
    const Epetra_CrsMatrix &mat, const TranspositionMode trans,
    const Eigen::Ref<Vector<std::complex<float>>> &x_in,
    Eigen::Ref<Vector<std::complex<float>>> y_inout, const std::complex<float> alpha,
    const std::complex<float> beta) {
  // Do the y_inout *= beta part
  const std::complex<float> zero(0.f, 0.f);
  if (beta == zero)
    y_inout.setZero();
  else
    y_inout *= beta;

  // Separate the real and imaginary components and store them in
  // double-precision vectors
  Vector<double> x_real(x_in.rows());
  for (size_t i = 0; i < x_in.rows(); ++i)
    x_real(i) = x_in(i).real();
  Vector<double> x_imag(x_in.rows());
  for (size_t i = 0; i < x_in.rows(); ++i)
    x_imag(i) = x_in(i).imag();
  Vector<double> y_real(y_inout.rows());
  for (size_t i = 0; i < y_inout.rows(); ++i)
    y_real(i) = y_inout(i).real();
  Vector<double> y_imag(y_inout.rows());
  for (size_t i = 0; i < y_inout.rows(); ++i)
    y_imag(i) = y_inout(i).imag();

  // Do the "+= alpha A x" part (in steps)
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real), Eigen::Ref<Vector<double>>(y_real), alpha.real(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag), Eigen::Ref<Vector<double>>(y_real), -alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real), Eigen::Ref<Vector<double>>(y_imag), alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag), Eigen::Ref<Vector<double>>(y_imag), alpha.real(), 1.);

  // Copy the result back to the complex vector
  for (size_t i = 0; i < y_inout.rows(); ++i)
    y_inout(i) = std::complex<float>(y_real(i), y_imag(i));
}

template <>
void reallyApplyBuiltInImpl<std::complex<double>>(
    const Epetra_CrsMatrix &mat, const TranspositionMode trans,
    const Eigen::Ref<Vector<std::complex<double>>> &x_in,
    Eigen::Ref<Vector<std::complex<double>>> y_inout, const std::complex<double> alpha,
    const std::complex<double> beta) {
  // Do the y_inout *= beta part
  const std::complex<double> zero(0., 0.);
  if (beta == zero)
    y_inout.fill(zero);
  else
    y_inout *= beta;

  // Separate the real and imaginary components
  Vector<double> x_real(x_in.real());
  Vector<double> x_imag(x_in.imag());
  Vector<double> y_real(y_inout.real());
  Vector<double> y_imag(y_inout.imag());

  // Do the "+= alpha A x" part (in steps)
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real), Eigen::Ref<Vector<double>>(y_real), alpha.real(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag), Eigen::Ref<Vector<double>>(y_real), -alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_real), Eigen::Ref<Vector<double>>(y_imag), alpha.imag(), 1.);
  reallyApplyBuiltInImpl<double>(mat, trans, Eigen::Ref<Vector<double>>(x_imag), Eigen::Ref<Vector<double>>(y_imag), alpha.real(), 1.);

  // Copy the result back to the complex vector
  for (size_t i = 0; i < y_inout.rows(); ++i)
    y_inout(i) = std::complex<double>(y_real(i), y_imag(i));
}

} // namespace

template <typename ValueType>
DiscreteSparseBoundaryOperator<ValueType>::DiscreteSparseBoundaryOperator(
    const shared_ptr<const Epetra_CrsMatrix> &mat, int symmetry,
    TranspositionMode trans, const shared_ptr<AhmedBemBlcluster> &blockCluster,
    const shared_ptr<IndexPermutation> &domainPermutation,
    const shared_ptr<IndexPermutation> &rangePermutation)
    : m_mat(mat), m_symmetry(symmetry), m_trans(trans),
      m_blockCluster(blockCluster), m_domainPermutation(domainPermutation),
      m_rangePermutation(rangePermutation) {
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::dump() const {
  if (isTransposed())
    std::cout << "Transpose of " << *m_mat << std::endl;
  else
    std::cout << *m_mat << std::endl;
}

template <typename ValueType>
Matrix<ValueType>
DiscreteSparseBoundaryOperator<ValueType>::asMatrix() const {
  if (m_mat->Comm().NumProc() != 1)
    throw std::runtime_error(
        "DiscreteSparseBoundaryOperator::asMatrix(): "
        "conversion of distributed matrices to local matrices is unsupported");

  bool transposed = isTransposed();
  const int untransposedRowCount = m_mat->NumGlobalRows();
  Matrix<ValueType> mat(rowCount(), columnCount());
  mat.setZero();
  for (int row = 0; row < untransposedRowCount; ++row) {
    int entryCount = 0;
    double *values = 0;
    int *indices = 0;
    int errorCode = m_mat->ExtractMyRowView(row, entryCount, values, indices);
    if (errorCode != 0)
      throw std::runtime_error("DiscreteSparseBoundaryOperator::asMatrix(): "
                               "Epetra_CrsMatrix::ExtractMyRowView()) failed");
    if (transposed)
      for (int entry = 0; entry < entryCount; ++entry)
        mat(indices[entry], row) = values[entry];
    else
      for (int entry = 0; entry < entryCount; ++entry)
        mat(row, indices[entry]) = values[entry];
  }
  return mat;
}

template <typename ValueType>
unsigned int DiscreteSparseBoundaryOperator<ValueType>::rowCount() const {
  return isTransposed() ? m_mat->NumGlobalCols() : m_mat->NumGlobalRows();
}

template <typename ValueType>
unsigned int DiscreteSparseBoundaryOperator<ValueType>::columnCount() const {
  return isTransposed() ? m_mat->NumGlobalRows() : m_mat->NumGlobalCols();
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  // indices of entries of the untransposed (stored) matrix
  bool transposed = isTransposed();
  const std::vector<int> &untransposedRows = transposed ? cols : rows;
  const std::vector<int> &untransposedCols = transposed ? rows : cols;

  if (block.rows() != rows.size() || block.cols() != cols.size())
    throw std::invalid_argument("DiscreteSparseBoundaryOperator::addBlock(): "
                                "incorrect block size");

  int entryCount = 0;
  double *values = 0;
  int *indices = 0;

  for (size_t row = 0; row < untransposedRows.size(); ++row) {
    // Provision for future MPI support.
    if (m_mat->IndicesAreLocal()) {
      int errorCode = m_mat->ExtractMyRowView(untransposedRows[row], entryCount,
                                              values, indices);
      if (errorCode != 0)
        throw std::runtime_error(
            "DiscreteSparseBoundaryOperator::addBlock(): "
            "Epetra_CrsMatrix::ExtractMyRowView()) failed");
    } else {
      int errorCode = m_mat->ExtractGlobalRowView(untransposedRows[row],
                                                  entryCount, values, indices);
      if (errorCode != 0)
        throw std::runtime_error(
            "DiscreteSparseBoundaryOperator::addBlock(): "
            "Epetra_CrsMatrix::ExtractGlobalRowView()) failed");
    }

    for (size_t col = 0; col < untransposedCols.size(); ++col)
      for (int entry = 0; entry < entryCount; ++entry)
        if (indices[entry] == cols[col])
          block(transposed ? col : row, transposed ? row : col) +=
              alpha * static_cast<ValueType>(values[entry]);
  }
}

#ifdef WITH_AHMED
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
DiscreteSparseBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
    double eps, int maximumRank, bool interleave) const {
  if (eps < 0)
    eps = 1e-4; // probably isn't really used

  // Note: the maximumRank parameter is ignored in this implementation.

  if (!m_blockCluster || !m_domainPermutation || !m_rangePermutation)
    throw std::runtime_error("DiscreteSparseBoundaryOperator::"
                             "asDiscreteAcaBoundaryOperator(): "
                             "data required for conversion to ACA operator "
                             "were not provided in the constructor");
  if (isTransposed())
    throw std::runtime_error("DiscreteSparseBoundaryOperator::"
                             "asDiscreteAcaBoundaryOperator(): "
                             "transposed operators are not supported yet");

  int *rowOffsets = 0;
  int *colIndices = 0;
  double *values = 0;
  m_mat->ExtractCrsDataPointers(rowOffsets, colIndices, values);

  std::vector<unsigned int> domain_o2p = m_domainPermutation->permutedIndices();
  std::vector<unsigned int> range_o2p = m_rangePermutation->permutedIndices();
  std::vector<unsigned int> range_p2o(range_o2p.size());
  for (size_t o = 0; o < range_o2p.size(); ++o) {
    assert(range_o2p[o] < range_o2p.size());
    range_p2o[range_o2p[o]] = o;
  }

  boost::shared_array<AhmedMblock *> mblocks;
  int trueMaximumRank = 0;
  SparseToHMatrixConverter<ValueType>::constructHMatrix(
      rowOffsets, colIndices, values, domain_o2p, range_p2o, eps,
      m_blockCluster.get(), mblocks, trueMaximumRank);

  // Gather remaining data necessary to create the combined ACA operator
  const int symmetry = 0;
  Fiber::ParallelizationOptions parallelOptions; // default options

  shared_ptr<const DiscreteBoundaryOperator<ValueType>> result(
      new DiscreteAcaBoundaryOperator<ValueType>(
          rowCount(), columnCount(), eps, trueMaximumRank, symmetry,
          m_blockCluster, mblocks, *m_domainPermutation, *m_rangePermutation,
          parallelOptions));
  return result;
}
#endif // WITH_AHMED

template <typename ValueType>
shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>>
DiscreteSparseBoundaryOperator<ValueType>::castToSparse(const shared_ptr<
    const DiscreteBoundaryOperator<ValueType>> &discreteOperator) {
  shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>> result =
      boost::dynamic_pointer_cast<
          const DiscreteSparseBoundaryOperator<ValueType>>(discreteOperator);
  if (result.get() == 0 && (discreteOperator.get() != 0))
    throw std::bad_cast();
  return result;
}

template <typename ValueType>
shared_ptr<const Epetra_CrsMatrix>
DiscreteSparseBoundaryOperator<ValueType>::epetraMatrix() const {
  return m_mat;
}

template <typename ValueType>
TranspositionMode
DiscreteSparseBoundaryOperator<ValueType>::transpositionMode() const {
  return m_trans;
}

template <typename ValueType>
bool DiscreteSparseBoundaryOperator<ValueType>::isTransposed() const {
  return m_trans & (TRANSPOSE | CONJUGATE_TRANSPOSE);
}

template <typename ValueType>
void DiscreteSparseBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {
  TranspositionMode realTrans = trans;
  bool transposed = isTransposed();
  if (transposed)
    switch (trans) {
    case NO_TRANSPOSE:
      realTrans = TRANSPOSE;
      break;
    case TRANSPOSE:
      realTrans = NO_TRANSPOSE;
      break;
    case CONJUGATE:
      realTrans = CONJUGATE_TRANSPOSE;
      break;
    case CONJUGATE_TRANSPOSE:
      realTrans = CONJUGATE;
      break;
      // default: should not happen; anyway, don't change trans
    }

  reallyApplyBuiltInImpl(*m_mat, realTrans, x_in, y_inout, alpha, beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteSparseBoundaryOperator);

} // namespace Bempp

