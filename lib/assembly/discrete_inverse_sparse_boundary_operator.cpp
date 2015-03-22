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
#ifdef WITH_TRILINOS

#include "discrete_inverse_sparse_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../common/eigen_support.hpp"

#include <iostream>
#include <stdexcept>

#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_SerialComm.h>
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>

namespace Bempp {

// Helper functions for the applyBuiltIn member function
namespace {

// Type-agnostic wrapper for the Amesos solver
template <typename ValueType>
void solveWithAmesos(Epetra_LinearProblem &problem, Amesos_BaseSolver &solver,
                     Matrix<ValueType> &solution,
                     const Matrix<ValueType> &rhs);

template <>
void solveWithAmesos<double>(Epetra_LinearProblem &problem,
                             Amesos_BaseSolver &solver,
                             Matrix<double> &armaSolution,
                             const Matrix<double> &armaRhs) {

  const size_t rowCount = armaRhs.rows();
  assert(rowCount == armaSolution.rows());
  const size_t rhsCount = armaRhs.cols();
  assert(rhsCount == armaSolution.cols());

  Epetra_Map map((int)rowCount, 0 /* base index */, Epetra_SerialComm());
  Epetra_MultiVector solution(View, map, armaSolution.memptr(), rowCount,
                              rhsCount);
  Epetra_MultiVector rhs(View, map, const_cast<double *>(armaRhs.memptr()),
                         rowCount, rhsCount);
  problem.SetLHS(&solution);
  problem.SetRHS(&rhs);

  if (solver.Solve() != 0)
    throw std::runtime_error("solveWithAmesos(): solve failed");
}

template <>
void solveWithAmesos<float>(Epetra_LinearProblem &problem,
                            Amesos_BaseSolver &solver,
                            Matrix<float> &armaSolution,
                            const Matrix<float> &armaRhs) {
  // Right now we only support single rhs vectors
  assert(armaSolution.cols() == 1);
  assert(armaRhs.cols() == 1);

  Vector<double> solution_double(armaSolution.rows());
  for (int i = 0; i < armaSolution.rows(); ++i) solution_double(i) = armaSolution(i,1);

  Vector<double> rhs_double(armaRhs.rows());
  for (int i = 0; i < armaRhs.rows(); ++i) rhs_double(i) = armaRhs(i,1);

  solveWithAmesos<double>(problem, solver, solution_double, rhs_double);

  for (int i = 0; i < solution_double.rows(); ++i) armaSolution(i,1) = solution_double(i);

}

template <>
void solveWithAmesos<std::complex<float>>(
    Epetra_LinearProblem &problem, Amesos_BaseSolver &solver,
    Matrix<std::complex<float>> &armaSolution,
    const Matrix<std::complex<float>> &armaRhs) {
  // Right now we only support single rhs vectors
  assert(armaSolution.cols() == 1);
  assert(armaRhs.cols() == 1);

  // Solve for the real and imaginary part separately
  // (The copy of the solution (before solving) is probably not necessary...)
  Matrix<double> solution_double(armaSolution.rows(), 2);
  for (size_t i = 0; i < armaSolution.rows(); ++i) {
    solution_double(i, 0) = armaSolution(i).real();
    solution_double(i, 1) = armaSolution(i).imag();
  }
  arma::Mat<double> rhs_double(armaRhs.rows(), 2);
  for (size_t i = 0; i < armaRhs.rows(); ++i) {
    rhs_double(i, 0) = armaRhs(i).real();
    rhs_double(i, 1) = armaRhs(i).imag();
  }

  solveWithAmesos<double>(problem, solver, solution_double, rhs_double);
  for (size_t i = 0; i < armaSolution.rows(); ++i)
    armaSolution(i) =
        std::complex<float>(solution_double(i, 0), solution_double(i, 1));
}

template <>
void solveWithAmesos<std::complex<double>>(
    Epetra_LinearProblem &problem, Amesos_BaseSolver &solver,
    Matrix<std::complex<double>> &armaSolution,
    const Matrix<std::complex<double>> &armaRhs) {
  // Right now we only support single rhs vectors
  assert(armaSolution.cols() == 1);
  assert(armaRhs.cols() == 1);

  // Solve for the real and imaginary part separately
  // (The copy of the solution (before solving) is probably not necessary...)
  arma::Mat<double> solution_double(armaSolution.rows(), 2);
  for (size_t i = 0; i < armaSolution.rows(); ++i) {
    solution_double(i, 0) = armaSolution(i).real();
    solution_double(i, 1) = armaSolution(i).imag();
  }
  arma::Mat<double> rhs_double(armaRhs.rows(), 2);
  for (size_t i = 0; i < armaRhs.rows(); ++i) {
    rhs_double(i, 0) = armaRhs(i).real();
    rhs_double(i, 1) = armaRhs(i).imag();
  }

  solveWithAmesos<double>(problem, solver, solution_double, rhs_double);
  for (size_t i = 0; i < armaSolution.rows(); ++i)
    armaSolution(i) =
        std::complex<double>(solution_double(i, 0), solution_double(i, 1));
}

} // namespace

template <typename ValueType>
DiscreteInverseSparseBoundaryOperator<ValueType>::
    DiscreteInverseSparseBoundaryOperator(
        const shared_ptr<const Epetra_CrsMatrix> &mat, int symmetry)
    : m_mat(mat), m_problem(new Epetra_LinearProblem),
      m_space(Thyra::defaultSpmdVectorSpace<ValueType>(mat->NumGlobalRows())),
      m_symmetry(symmetry) {

  if (m_mat->NumGlobalRows() != m_mat->NumGlobalCols())
    throw std::invalid_argument("DiscreteInverseSparseBoundaryOperator::"
                                "DiscreteInverseSparseBoundaryOperator(): "
                                "square matrix expected");
  // const_cast: Amesos is not const-correct. Amesos2 will be,
  // and Amesos2 takes a RCP to a const matrix.
  m_problem->SetOperator(const_cast<Epetra_CrsMatrix *>(m_mat.get()));
  if (m_symmetry & (SYMMETRIC | HERMITIAN)) // Epetra matrices are real, so
                                            // symmetric == Hermitian
    m_problem->AssertSymmetric();

  Amesos amesosFactory;
  const char *solverName = "Amesos_Klu";
  if (!amesosFactory.Query(solverName))
    throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                             "DiscreteInverseSparseBoundaryOperator(): "
                             "Amesos_Klu solver not available");
  m_solver.reset(amesosFactory.Create("Amesos_Klu", *m_problem));

  if (!m_solver.get())
    throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                             "DiscreteInverseSparseBoundaryOperator(): "
                             "Amesos solver could not be constructed");
  if (m_solver->SymbolicFactorization() != 0)
    throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                             "DiscreteInverseSparseBoundaryOperator(): "
                             "Symbolic factorization with Amesos failed");
  if (m_solver->NumericFactorization() != 0)
    throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                             "DiscreteInverseSparseBoundaryOperator(): "
                             "Numeric factorization with Amesos failed");
}

template <typename ValueType>
DiscreteInverseSparseBoundaryOperator<
    ValueType>::~DiscreteInverseSparseBoundaryOperator() {}

template <typename ValueType>
unsigned int
DiscreteInverseSparseBoundaryOperator<ValueType>::rowCount() const {
  return m_space->dim();
}

template <typename ValueType>
unsigned int
DiscreteInverseSparseBoundaryOperator<ValueType>::columnCount() const {
  return m_space->dim();
}

template <typename ValueType>
void DiscreteInverseSparseBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  throw std::runtime_error("DiscreteInverseSparseBoundaryOperator::"
                           "addBlock(): not implemented");
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteInverseSparseBoundaryOperator<ValueType>::domain() const {
  return m_space;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteInverseSparseBoundaryOperator<ValueType>::range() const {
  return m_space;
}

template <typename ValueType>
bool DiscreteInverseSparseBoundaryOperator<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  return (M_trans == Thyra::NOTRANS);
}

template <typename ValueType>
void DiscreteInverseSparseBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  // TODO: protect with a mutex (this function is not thread-safe)
  if (trans != NO_TRANSPOSE)
    throw std::invalid_argument("DiscreteInverseSparseBoundaryOperator::"
                                "applyBuiltInImpl(): "
                                "transposes and conjugates are not supported");
  const size_t dim = m_space->dim();
  if (x_in.rows() != dim || y_inout.rows() != dim)
    throw std::invalid_argument("DiscreteInverseSparseBoundaryOperator::"
                                "applyBuiltInImpl(): "
                                "incorrect vector lengths");
  Vector<ValueType> solution(dim);
  solution.setZero();
  solveWithAmesos(*m_problem, *m_solver, solution, x_in);
  if (beta == static_cast<ValueType>(0.))
    y_inout = alpha * solution;
  else {
    y_inout *= beta;
    y_inout += alpha * solution;
  }
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>> discreteSparseInverse(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &discreteOp) {
  shared_ptr<const DiscreteSparseBoundaryOperator<ValueType>> sparseOp =
      DiscreteSparseBoundaryOperator<ValueType>::castToSparse(discreteOp);

  shared_ptr<const DiscreteBoundaryOperator<ValueType>> op(
      new DiscreteInverseSparseBoundaryOperator<ValueType>(
          sparseOp->epetraMatrix(), sparseOp->symmetryMode()));
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

#endif // WITH_TRILINOS
