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

#include "boundary_operator.hpp"

#include "abstract_boundary_operator_sum.hpp"
#include "abstract_boundary_operator_composition.hpp"
#include "adjoint_abstract_boundary_operator.hpp"
#include "discrete_boundary_operator.hpp"
#include "context.hpp"
#include "grid_function.hpp"
#include "scaled_abstract_boundary_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>::BoundaryOperator()
    : m_holdWeakForm(true) {}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>::~BoundaryOperator() {}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>::BoundaryOperator(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const AbstractBoundaryOperator<BasisFunctionType,
                                                    ResultType>> &abstractOp)
    : m_holdWeakForm(true) {
  initialize(context, abstractOp);
}

template <typename BasisFunctionType, typename ResultType>
void BoundaryOperator<BasisFunctionType, ResultType>::initialize(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const AbstractBoundaryOperator<BasisFunctionType,
                                                    ResultType>> &abstractOp) {
  if (!context)
    throw std::invalid_argument("BoundaryOperator::BoundaryOperator(): "
                                "context must not be null");
  if (!abstractOp)
    throw std::invalid_argument("BoundaryOperator::BoundaryOperator(): "
                                "abstractOp must not be null");
  m_context = context;
  m_abstractOp = abstractOp;
  m_weakFormContainer.reset(new ConstWeakFormContainer);
  m_weakWeakFormContainer.reset(new WeakConstWeakFormContainer);
}

template <typename BasisFunctionType, typename ResultType>
void BoundaryOperator<BasisFunctionType, ResultType>::uninitialize() {
  m_context.reset();
  m_abstractOp.reset();
  m_weakFormContainer.reset();
  m_weakWeakFormContainer.reset();
}

template <typename BasisFunctionType, typename ResultType>
bool BoundaryOperator<BasisFunctionType, ResultType>::isInitialized() const {
  return (bool)m_abstractOp;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType>>
BoundaryOperator<BasisFunctionType, ResultType>::abstractOperator() const {
  return m_abstractOp;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Context<BasisFunctionType, ResultType>>
BoundaryOperator<BasisFunctionType, ResultType>::context() const {
  return m_context;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
BoundaryOperator<BasisFunctionType, ResultType>::weakForm(bool update) const {
  if (!isInitialized())
    throw std::runtime_error(
        "BoundaryOperator::weakForm(): attempted to retrieve the "
        "weak form of an uninitialized operator");
  assert(m_weakFormContainer);     // contains a shared_ptr to DiscreteOp
                                   // (which may be null, though)
  assert(m_weakWeakFormContainer); // contains a weak_ptr to DiscreteOp
                                   // (which may be null, though)
  typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
  shared_ptr<const DiscreteOp> discreteOp = m_weakWeakFormContainer->lock();
  if (!discreteOp || update) {
    discreteOp = m_abstractOp->assembleWeakForm(*m_context);
    assert(discreteOp);
    *m_weakWeakFormContainer = discreteOp;
    if (m_holdWeakForm)
      *m_weakFormContainer = discreteOp;
  }
  return discreteOp;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
BoundaryOperator<BasisFunctionType, ResultType>::domain() const {
  if (!isInitialized())
    return shared_ptr<const Space<BasisFunctionType>>();
  return m_abstractOp->domain();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
BoundaryOperator<BasisFunctionType, ResultType>::range() const {
  if (!isInitialized())
    return shared_ptr<const Space<BasisFunctionType>>();
  return m_abstractOp->range();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
BoundaryOperator<BasisFunctionType, ResultType>::dualToRange() const {
  if (!isInitialized())
    return shared_ptr<const Space<BasisFunctionType>>();
  return m_abstractOp->dualToRange();
}

template <typename BasisFunctionType, typename ResultType>
std::string BoundaryOperator<BasisFunctionType, ResultType>::label() const {
  if (!isInitialized())
    return std::string();
  return m_abstractOp->label();
}

template <typename BasisFunctionType, typename ResultType>
bool BoundaryOperator<BasisFunctionType, ResultType>::isWeakFormHeld() const {
  return m_holdWeakForm;
}

template <typename BasisFunctionType, typename ResultType>
void BoundaryOperator<BasisFunctionType, ResultType>::holdWeakForm(bool value) {
  if (value)
    *m_weakFormContainer = m_weakWeakFormContainer->lock();
  else {
    m_weakFormContainer->reset();
  }
  m_holdWeakForm = value;
}

template <typename BasisFunctionType, typename ResultType>
void BoundaryOperator<BasisFunctionType, ResultType>::apply(
    const TranspositionMode trans,
    const GridFunction<BasisFunctionType, ResultType> &x_in,
    GridFunction<BasisFunctionType, ResultType> &y_inout, ResultType alpha,
    ResultType beta) const {
  if (!isInitialized())
    throw std::runtime_error("BoundaryOperator::apply(): attempted to apply "
                             "an uninitialized operator");

  // Sanity test
  if (!m_abstractOp->domain()->spaceIsCompatible(*x_in.space()) ||
      !m_abstractOp->range()->spaceIsCompatible(*y_inout.space()))
    throw std::invalid_argument("BoundaryOperator::apply(): "
                                "spaces don't match");

  shared_ptr<const Space<BasisFunctionType>> dualToRange =
      m_abstractOp->dualToRange();

  // Extract coefficient vectors
  Matrix<ResultType> xVals = x_in.coefficients();
  Matrix<ResultType> yVals = y_inout.projections(dualToRange);

  // Apply operator and assign the result to y_inout's projections
  weakForm()->apply(trans, xVals, yVals, alpha, beta);
  // TODO: make interfaces to the Trilinos and fallback
  // DiscreteBoundaryOperator::apply() compatible.
  // Perhaps by declaring an asPtrToBaseVector method in Vector...
  y_inout.setProjections(dualToRange, yVals.col(0));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator+(const BoundaryOperator<BasisFunctionType, ResultType> &op) {
  return op;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator-(const BoundaryOperator<BasisFunctionType, ResultType> &op) {
  if (!op.isInitialized())
    throw std::invalid_argument("operator-(): operand is uninitialized");
  return static_cast<ResultType>(-1.) * op;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator+(const BoundaryOperator<BasisFunctionType, ResultType> &op1,
          const BoundaryOperator<BasisFunctionType, ResultType> &op2) {
  typedef AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> Sum;
  if (!op1.isInitialized())
    throw std::invalid_argument("operator+(): operand 1 is uninitialized");
  if (!op2.isInitialized())
    throw std::invalid_argument("operator+(): operand 2 is uninitialized");
  return BoundaryOperator<BasisFunctionType, ResultType>(
      op1.context(), boost::make_shared<Sum>(op1, op2));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator-(const BoundaryOperator<BasisFunctionType, ResultType> &op1,
          const BoundaryOperator<BasisFunctionType, ResultType> &op2) {
  if (!op1.isInitialized())
    throw std::invalid_argument("operator-(): operand 1 is uninitialized");
  if (!op2.isInitialized())
    throw std::invalid_argument("operator-(): operand 2 is uninitialized");
  return op1 + (static_cast<ResultType>(-1.) * op2);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>,
                        std::complex<double>>,
        ScalarType>,
    BoundaryOperator<BasisFunctionType, ResultType>>::type
operator*(const BoundaryOperator<BasisFunctionType, ResultType> &op,
          const ScalarType &scalar) {
  typedef ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType>
      ScaledOp;
  if (!op.isInitialized())
    throw std::invalid_argument("operator*(): "
                                "boundary operator is uninitialized");
  return BoundaryOperator<BasisFunctionType, ResultType>(
      op.context(),
      boost::make_shared<ScaledOp>(static_cast<ResultType>(scalar), op));
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType>
operator*(const ScalarType &scalar,
          const BoundaryOperator<BasisFunctionType, ResultType> &op) {
  return operator*(op, scalar);
}

template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType>
operator/(const BoundaryOperator<BasisFunctionType, ResultType> &op,
          const ScalarType &scalar) {
  if (!op.isInitialized())
    throw std::runtime_error("operator/(BoundaryOperator, scalar): "
                             "operand is uninitialized");
  if (scalar == static_cast<ScalarType>(0.))
    throw std::runtime_error("operator/(BoundaryOperator, scalar): "
                             "Division by zero");
  return operator*(
      op, static_cast<ResultType>(static_cast<ScalarType>(1.) / scalar));
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
operator*(const BoundaryOperator<BasisFunctionType, ResultType> &op,
          const GridFunction<BasisFunctionType, ResultType> &fun) {
  if (!op.isInitialized())
    throw std::runtime_error("operator*(BoundaryOperator, GridFunction): "
                             "operand 1 is uninitialized");
  if (!fun.isInitialized())
    throw std::runtime_error("operator*(BoundaryOperator, GridFunction): "
                             "operand 2 is uninitialized");

  typedef GridFunction<BasisFunctionType, ResultType> GF;

  shared_ptr<const Space<BasisFunctionType>> space = op.range();
  shared_ptr<const Space<BasisFunctionType>> dualSpace = op.dualToRange();
  Vector<ResultType> projections(dualSpace->globalDofCount());
  projections.setZero();
  GF result(op.context(), space, dualSpace, projections);
  op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
  return result;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator*(const BoundaryOperator<BasisFunctionType, ResultType> &op1,
          const BoundaryOperator<BasisFunctionType, ResultType> &op2) {
  typedef AbstractBoundaryOperatorComposition<BasisFunctionType, ResultType>
      Composition;
  if (!op1.isInitialized())
    throw std::invalid_argument("operator*(): operand 1 is uninitialized");
  if (!op2.isInitialized())
    throw std::invalid_argument("operator*(): operand 2 is uninitialized");
  return BoundaryOperator<BasisFunctionType, ResultType>(
      op1.context(), boost::make_shared<Composition>(op1, op2));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
adjoint(const BoundaryOperator<BasisFunctionType, ResultType> &op) {
  typedef AdjointAbstractBoundaryOperator<BasisFunctionType, ResultType>
      Adjoint;
  return BoundaryOperator<BasisFunctionType, ResultType>(
      op.context(), boost::make_shared<Adjoint>(op));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
adjoint(const BoundaryOperator<BasisFunctionType, ResultType> &op,
        const shared_ptr<const Space<BasisFunctionType>> &range) {
  typedef AdjointAbstractBoundaryOperator<BasisFunctionType, ResultType>
      Adjoint;
  return BoundaryOperator<BasisFunctionType, ResultType>(
      op.context(), boost::make_shared<Adjoint>(op, range));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> &
throwIfUninitialized(BoundaryOperator<BasisFunctionType, ResultType> &op,
                     std::string message) {
  if (op.isInitialized())
    return op;
  else if (message.empty())
    message = "Detected an uninitialized BoundaryOperator object";
  throw std::invalid_argument(message);
}

template <typename BasisFunctionType, typename ResultType>
const BoundaryOperator<BasisFunctionType, ResultType> &
throwIfUninitialized(const BoundaryOperator<BasisFunctionType, ResultType> &op,
                     std::string message) {
  if (op.isInitialized())
    return op;
  else if (message.empty())
    message = "Detected an uninitialized BoundaryOperator object";
  throw std::invalid_argument(message);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(BASIS, RESULT)                              \
  template BoundaryOperator<BASIS, RESULT> operator+(                          \
      const BoundaryOperator<BASIS, RESULT> &op);                              \
  template BoundaryOperator<BASIS, RESULT> operator-(                          \
      const BoundaryOperator<BASIS, RESULT> &op);                              \
  template BoundaryOperator<BASIS, RESULT> operator+(                          \
      const BoundaryOperator<BASIS, RESULT> &op1,                              \
      const BoundaryOperator<BASIS, RESULT> &op2);                             \
  template BoundaryOperator<BASIS, RESULT> operator-(                          \
      const BoundaryOperator<BASIS, RESULT> &op1,                              \
      const BoundaryOperator<BASIS, RESULT> &op2);                             \
  template GridFunction<BASIS, RESULT> operator*(                              \
      const BoundaryOperator<BASIS, RESULT> &op,                               \
      const GridFunction<BASIS, RESULT> &fun);                                 \
  template BoundaryOperator<BASIS, RESULT> operator*(                          \
      const BoundaryOperator<BASIS, RESULT> &op,                               \
      const BoundaryOperator<BASIS, RESULT> &fun);                             \
  template BoundaryOperator<BASIS, RESULT> adjoint(                            \
      const BoundaryOperator<BASIS, RESULT> &op);                              \
  template BoundaryOperator<BASIS, RESULT> adjoint(                            \
      const BoundaryOperator<BASIS, RESULT> &op,                               \
      const shared_ptr<const Space<BASIS>> &range);                            \
  template BoundaryOperator<BASIS, RESULT> &throwIfUninitialized(              \
      BoundaryOperator<BASIS, RESULT> &op, std::string message);               \
  template const BoundaryOperator<BASIS, RESULT> &throwIfUninitialized(        \
      const BoundaryOperator<BASIS, RESULT> &op, std::string message)
#define INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(BASIS, RESULT, SCALAR)          \
  template BoundaryOperator<BASIS, RESULT> operator*(                          \
      const BoundaryOperator<BASIS, RESULT> &op, const SCALAR &scalar);        \
  template BoundaryOperator<BASIS, RESULT> operator*(                          \
      const SCALAR &scalar, const BoundaryOperator<BASIS, RESULT> &op);        \
  template BoundaryOperator<BASIS, RESULT> operator/(                          \
      const BoundaryOperator<BASIS, RESULT> &op, const SCALAR &scalar)

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, float, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, float, double);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
INSTANTIATE_FREE_FUNCTIONS(float, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, std::complex<float>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, std::complex<float>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, std::complex<float>,
                                       std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(float, std::complex<float>,
                                       std::complex<double>);
#endif

#if defined(ENABLE_SINGLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
INSTANTIATE_FREE_FUNCTIONS(std::complex<float>, std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, std::complex<float>,
                                       float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, std::complex<float>,
                                       double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, std::complex<float>,
                                       std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<float>, std::complex<float>,
                                       std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(double, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, double, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, double, double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_KERNELS)
INSTANTIATE_FREE_FUNCTIONS(double, std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, std::complex<double>,
                                       std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(double, std::complex<double>,
                                       std::complex<double>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) && defined(ENABLE_COMPLEX_BASIS_FUNCTIONS)
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>, std::complex<double>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>,
                                       std::complex<double>, float);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>,
                                       std::complex<double>, double);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>,
                                       std::complex<double>,
                                       std::complex<float>);
INSTANTIATE_FREE_FUNCTIONS_WITH_SCALAR(std::complex<double>,
                                       std::complex<double>,
                                       std::complex<double>);
#endif

} // namespace Bempp
