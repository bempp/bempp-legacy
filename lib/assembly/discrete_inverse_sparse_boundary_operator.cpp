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


#include "discrete_inverse_sparse_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../common/eigen_support.hpp"

#include <iostream>
#include <stdexcept>

namespace {

template <typename ValueType>
Bempp::Vector<ValueType> solveWithEigen(const Eigen::SparseLU<Bempp::RealSparseMatrix>& solver,
                                 const Eigen::Ref<Bempp::Vector<ValueType>>& x);

template <>
Bempp::Vector<double> solveWithEigen<double>(const Eigen::SparseLU<Bempp::RealSparseMatrix>& solver,
                                 const Eigen::Ref<Bempp::Vector<double>>& x){

    return solver.solve(x);

}

template <>
Bempp::Vector<float> solveWithEigen<float>(const Eigen::SparseLU<Bempp::RealSparseMatrix>& solver,
                                 const Eigen::Ref<Bempp::Vector<float>>& x){

    return solver.solve(x.template cast<double>()).template cast<float>();

}

template <>
Bempp::Vector<std::complex<float>> solveWithEigen<std::complex<float>>(const Eigen::SparseLU<Bempp::RealSparseMatrix>& solver,
                                 const Eigen::Ref<Bempp::Vector<std::complex<float>>>& x){

    Bempp::Vector<double> x_re = x.real().template cast<double>();
    Bempp::Vector<double> x_im = x.imag().template cast<double>();

    Bempp::Vector<std::complex<float>> result(x.rows());
    result.real() = solver.solve(x_re).template cast<float>();
    result.imag() = solver.solve(x_im).template cast<float>();
    return result;
}

template <>
Bempp::Vector<std::complex<double>> solveWithEigen<std::complex<double>>(const Eigen::SparseLU<Bempp::RealSparseMatrix>& solver,
                                 const Eigen::Ref<Bempp::Vector<std::complex<double>>>& x){

    Bempp::Vector<std::complex<double>> result(x.rows());
    Bempp::Vector<double> x_real = x.real();
    Bempp::Vector<double> x_imag = x.imag();

    Bempp::Vector<double> res_real = solver.solve(x_real);
    Bempp::Vector<double> res_imag = solver.solve(x_imag);

    result.real() = res_real;
    result.imag() = res_imag;
    return result;
}



}

namespace Bempp {



template <typename ValueType>
DiscreteInverseSparseBoundaryOperator<ValueType>::
    DiscreteInverseSparseBoundaryOperator(
        const shared_ptr<const RealSparseMatrix> &mat, int symmetry)
    : m_mat(mat), m_solver(new Eigen::SparseLU<RealSparseMatrix>()),
      m_symmetry(symmetry) {

  if (m_mat->rows() != m_mat->cols())
    throw std::invalid_argument("DiscreteInverseSparseBoundaryOperator::"
                                "DiscreteInverseSparseBoundaryOperator(): "
                                "square matrix expected");

  m_solver->compute(*mat);

  if (m_solver->info()!=Eigen::Success)
      throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                               "DiscreteInverseSparseBoundaryOperator(): "
                               "LU Decomposition failed.");

}

template <typename ValueType>
DiscreteInverseSparseBoundaryOperator<
    ValueType>::~DiscreteInverseSparseBoundaryOperator() {}

template <typename ValueType>
unsigned int
DiscreteInverseSparseBoundaryOperator<ValueType>::rowCount() const {
  return m_mat->rows();
}

template <typename ValueType>
unsigned int
DiscreteInverseSparseBoundaryOperator<ValueType>::columnCount() const {
  return m_mat->cols();
}

template <typename ValueType>
void DiscreteInverseSparseBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                           "addBlock(): not implemented");
}

template <typename ValueType>
void DiscreteInverseSparseBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {
  // TODO: protect with a mutex (this function is not thread-safe)
  if (trans != NO_TRANSPOSE)
    throw std::invalid_argument("DiscreteInverseSparseBoundaryOperator::"
                                "applyBuiltInImpl(): "
                                "transposes and conjugates are not supported");
  const size_t dim = m_mat->rows();
  if (x_in.rows() != dim || y_inout.rows() != dim)
    throw std::invalid_argument("DiscreteInverseSparseBoundaryOperator::"
                                "applyBuiltInImpl(): "
                                "incorrect vector lengths");

  if (beta == static_cast<ValueType>(0.)){
    y_inout = alpha * solveWithEigen<ValueType>(*m_solver,x_in);
  }
  else {
    y_inout *= beta;
    y_inout += alpha * solveWithEigen<ValueType>(*m_solver,x_in);
  }
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>> discreteSparseInverse(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &discreteOp) {
  shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>> sparseOp =
      DiscreteSparseBoundaryOperator<ValueType>::castToSparse(discreteOp);

  shared_ptr<const DiscreteBoundaryOperator<ValueType>> op(
      new DiscreteInverseSparseBoundaryOperator<ValueType>(
          sparseOp->sparseMatrix(), sparseOp->symmetryMode()));
  return op;
}

#define INSTANTIATE_FREE_FUNCTIONS(VALUE)                                      \
  template shared_ptr<const DiscreteBoundaryOperator<VALUE>>                   \
  discreteSparseInverse(                                                       \
      const shared_ptr<const DiscreteBoundaryOperator<VALUE>> &discreteOp);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(
    DiscreteInverseSparseBoundaryOperator);
FIBER_ITERATE_OVER_VALUE_TYPES(INSTANTIATE_FREE_FUNCTIONS);

} // namespace Bempp

