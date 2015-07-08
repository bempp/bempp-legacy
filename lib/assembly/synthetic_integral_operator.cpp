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

#include "synthetic_integral_operator.hpp"
#include "../common/eigen_support.hpp"

#include "boundary_operator.hpp"
#include "context.hpp"
#include "discrete_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"
#include "identity_operator.hpp"
#include "sparse_inverse.hpp"
#include "transposed_discrete_boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <boost/make_shared.hpp>
#include <tbb/tick_count.h>

namespace Bempp {

namespace {

template <typename T> const T &vectorFirstElement(const std::vector<T> &v) {
  if (v.empty())
    throw std::invalid_argument(
        "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
        "lists of test and trial-side operators must not be empty");
  return v.front();
}

template <typename ResultType>
std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
coalesceTestOperators(
    const std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
        &discreteLocalOps,
    const shared_ptr<const RealSparseMatrix> &idInverse) {
  typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
  typedef DiscreteSparseBoundaryOperator<ResultType> SparseOp;
  std::vector<shared_ptr<const DiscreteOp>> result(discreteLocalOps.size());

  for (size_t i = 0; i < discreteLocalOps.size(); ++i) {
    if (!discreteLocalOps[i])
      throw std::runtime_error(
          "coalesceTrialOperators(): null pointer to operator detected");
    shared_ptr<const SparseOp> op =
        boost::dynamic_pointer_cast<const SparseOp>(discreteLocalOps[i]);
    if (op) {
      shared_ptr<const RealSparseMatrix> opMat = op->sparseMatrix();
      shared_ptr<RealSparseMatrix> composition(
          new RealSparseMatrix(*opMat * (*idInverse)));
      result[i].reset(new SparseOp(composition));
    } else
      throw std::runtime_error(
          "SyntheticIntegralOperator::coalesceTestOperators(): "
          "local operator not represented by a sparse matrix");
  }
  return result;
}

template <typename ResultType>
std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
coalesceTrialOperators(
    const std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
        &discreteLocalOps,
    const shared_ptr<const RealSparseMatrix> &idInverse) {
  typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
  typedef DiscreteSparseBoundaryOperator<ResultType> SparseOp;
  std::vector<shared_ptr<const DiscreteOp>> result(discreteLocalOps.size());

  for (size_t i = 0; i < discreteLocalOps.size(); ++i) {
    if (!discreteLocalOps[i])
      throw std::runtime_error(
          "coalesceTrialOperators(): null pointer to operator detected");
    shared_ptr<const SparseOp> op =
        boost::dynamic_pointer_cast<const SparseOp>(discreteLocalOps[i]);
    if (op) {
      shared_ptr<const RealSparseMatrix> opMat = op->sparseMatrix();
      shared_ptr<RealSparseMatrix> composition(
          new RealSparseMatrix(*idInverse * (*opMat)));
      result[i].reset(new SparseOp(composition));
    } else
      throw std::runtime_error(
          "SyntheticIntegralOperator::coalesceTrialOperators(): "
          "local operator not represented by a sparse matrix");
  }
  return result;
}

template <typename ResultType>
std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
transposeTestOperators(
    const std::vector<shared_ptr<const DiscreteBoundaryOperator<ResultType>>>
        &discreteLocalOps,
    bool hermitian) {
  typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
  typedef DiscreteSparseBoundaryOperator<ResultType> SparseOp;
  std::vector<shared_ptr<const DiscreteOp>> result(discreteLocalOps.size());

  for (size_t i = 0; i < discreteLocalOps.size(); ++i) {
    if (!discreteLocalOps[i])
      throw std::runtime_error(
          "coalesceTrialOperators(): null pointer to operator detected");
    shared_ptr<const SparseOp> op =
        boost::dynamic_pointer_cast<const SparseOp>(discreteLocalOps[i]);
    if (op) {
      int mode = op->transpositionMode();
      if (mode & TRANSPOSE)
        mode &= ~TRANSPOSE;
      else
        mode |= TRANSPOSE;
      result[i].reset(new SparseOp(op->sparseMatrix(), op->symmetryMode(),
                                   static_cast<TranspositionMode>(mode)));
    } else
      result[i].reset(new TransposedDiscreteBoundaryOperator<ResultType>(
          hermitian ? CONJUGATE_TRANSPOSE : TRANSPOSE, discreteLocalOps[i]));
  }
  return result;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>> determineDomain(
    const std::vector<BoundaryOperator<BasisFunctionType, ResultType>>
        &testLocalOps,
    const BoundaryOperator<BasisFunctionType, ResultType> &integralOp,
    const std::vector<BoundaryOperator<BasisFunctionType, ResultType>>
        &trialLocalOps,
    int symmetry) {
  if (trialLocalOps.empty())
    if (testLocalOps.empty())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
          "testLocalOps and trialLocalOps must not both be empty.");
    else {
      if (symmetry & (SYMMETRIC | HERMITIAN))
        return testLocalOps[0].dualToRange();
      else
        return integralOp.domain();
    }
  else
    return trialLocalOps[0].domain();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>> determineRange(
    const std::vector<BoundaryOperator<BasisFunctionType, ResultType>>
        &testLocalOps,
    const BoundaryOperator<BasisFunctionType, ResultType> &integralOp) {
  if (testLocalOps.empty())
    return integralOp.range();
  else
    return testLocalOps[0].range();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>> determineDualToRange(
    const std::vector<BoundaryOperator<BasisFunctionType, ResultType>>
        &testLocalOps,
    const BoundaryOperator<BasisFunctionType, ResultType> &integralOp) {
  if (testLocalOps.empty())
    return integralOp.dualToRange();
  else
    return testLocalOps[0].dualToRange();
}

} // namespace

template <typename BasisFunctionType, typename ResultType>
SyntheticIntegralOperator<BasisFunctionType, ResultType>::
    SyntheticIntegralOperator(
        const std::vector<BoundaryOperator<BasisFunctionType, ResultType>>
            &testLocalOps,
        const BoundaryOperator<BasisFunctionType, ResultType> &integralOp,
        const std::vector<BoundaryOperator<BasisFunctionType, ResultType>>
            &trialLocalOps,
        const std::string &label, int syntheseSymmetry)
    : Base(determineDomain(testLocalOps, integralOp, trialLocalOps,
                           syntheseSymmetry),
           determineRange(testLocalOps, integralOp),
           determineDualToRange(testLocalOps, integralOp), label,
           syntheseSymmetry & integralOp.abstractOperator()->symmetry()),
      m_integralOp(integralOp), m_testLocalOps(testLocalOps),
      m_trialLocalOps(trialLocalOps), m_syntheseSymmetry(syntheseSymmetry) {
  // Note: the code does not at present distinguish properly between symmetric
  // and/or hermitian operators (in fact symmetry handling is, frankly, a
  // mess). We get away with this because sparse operators can only contain
  // real entries anyway. To be fixed in future.

  checkIntegralOperator();
  if (m_testLocalOps.size() > 0 && m_trialLocalOps.size() > 0 &&
      m_testLocalOps.size() != m_trialLocalOps.size())
    throw std::invalid_argument(
        "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
        "if both testLocalOps and trialLocalOps are non-empty, "
        "both must have the same number of elements");
  checkTestLocalOperators();
  if ((syntheseSymmetry & SYMMETRIC) || (syntheseSymmetry & HERMITIAN)) {
    if (m_testLocalOps.empty())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
          "symmetric operator requested but testLocalOps empty");
    m_trialLocalOps.clear();
  } else
    checkTrialLocalOperators();
  if (m_testLocalOps.empty() && m_trialLocalOps.empty())
    throw std::invalid_argument(
        "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
        "testLocalOps and trialLocalOps must not both be empty.");
}

template <typename BasisFunctionType, typename ResultType>
bool SyntheticIntegralOperator<BasisFunctionType, ResultType>::isLocal() const {
  return false;
}

template <typename BasisFunctionType, typename ResultType>
void SyntheticIntegralOperator<BasisFunctionType,
                               ResultType>::checkIntegralOperator() const {
  if (!m_integralOp.isInitialized())
    throw std::invalid_argument(
        "SyntheticIntegralOperator::checkIntegralOperator(): "
        "all operators must be initialized");
}

template <typename BasisFunctionType, typename ResultType>
void SyntheticIntegralOperator<BasisFunctionType,
                               ResultType>::checkTestLocalOperators() const {
  size_t localOperatorCount = m_testLocalOps.size();
  for (size_t i = 0; i < localOperatorCount; ++i)
    if (!m_testLocalOps[i].isInitialized())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTestLocalOperators(): "
          "all operators must be initialized");
  for (size_t i = 0; i < localOperatorCount; ++i)
    if (!m_testLocalOps[i].abstractOperator()->isLocal())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTestLocalOperators(): "
          "all test operators must be local");
  for (size_t i = 0; i < localOperatorCount; ++i)
    if (m_testLocalOps[i].domain() != m_integralOp.dualToRange())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTestLocalOperators(): "
          "domain of all test operators must be identical with the dual "
          "to range of the integral operator");
  for (size_t i = 1; i < localOperatorCount; ++i)
    if (m_testLocalOps[i].domain() != m_testLocalOps[0].domain() ||
        m_testLocalOps[i].range() != m_testLocalOps[0].range() ||
        m_testLocalOps[i].dualToRange() != m_testLocalOps[0].dualToRange())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTestLocalOperators(): "
          "test operators have inconsistent spaces");
}

template <typename BasisFunctionType, typename ResultType>
void SyntheticIntegralOperator<BasisFunctionType,
                               ResultType>::checkTrialLocalOperators() const {
  size_t localOperatorCount = m_trialLocalOps.size();
  for (size_t i = 0; i < localOperatorCount; ++i)
    if (!m_trialLocalOps[i].isInitialized())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTrialLocalOperators(): "
          "all operators must be initialized");
  for (size_t i = 0; i < localOperatorCount; ++i)
    if (!m_trialLocalOps[i].abstractOperator()->isLocal())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTrialLocalOperators(): "
          "all test operators must be local");
  for (size_t i = 0; i < localOperatorCount; ++i)
    if (m_trialLocalOps[i].dualToRange() != m_integralOp.domain())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTrialLocalOperators(): "
          "dual to range of all trial operators must be identical with the "
          "domain of the integral operator");
  for (size_t i = 1; i < localOperatorCount; ++i)
    if (m_trialLocalOps[i].domain() != m_trialLocalOps[0].domain() ||
        m_trialLocalOps[i].range() != m_trialLocalOps[0].range() ||
        m_trialLocalOps[i].dualToRange() != m_trialLocalOps[0].dualToRange())
      throw std::invalid_argument(
          "SyntheticIntegralOperator::checkTrialLocalOperators(): "
          "trial operators have inconsistent spaces");
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
SyntheticIntegralOperator<BasisFunctionType, ResultType>::assembleWeakFormImpl(
    const Context<BasisFunctionType, ResultType> &context) const {
  bool verbose =
      (context.assemblyOptions().verbosityLevel() >= VerbosityLevel::DEFAULT);
  if (verbose)
    std::cout << "Assembling the weak form of operator '" << this->label()
              << "'..." << std::endl;

  tbb::tick_count start = tbb::tick_count::now();

  typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
  typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;

  bool symmetricMode =
      m_syntheseSymmetry & SYMMETRIC || m_syntheseSymmetry & HERMITIAN;
  // std::cout << "symmetricMode: " << symmetricMode << std::endl;

  // We don't need a persistent shared pointer
  shared_ptr<const Context<BasisFunctionType, ResultType>> internalContext,
      auxContext;
  getContextsForInternalAndAuxiliaryOperators(make_shared_from_ref(context),
                                              internalContext, auxContext);

  shared_ptr<const DiscreteLinOp> discreteIntegralOp = m_integralOp.weakForm();
  std::vector<shared_ptr<const DiscreteLinOp>> discreteTestLocalOps(
      m_testLocalOps.size()),
      discreteTrialLocalOps(m_trialLocalOps.size());
  for (size_t i = 0; i < m_testLocalOps.size(); ++i)
    discreteTestLocalOps[i] = m_testLocalOps[i].weakForm();
  // doesn't do any harm if the list of trial operators is empty
  for (size_t i = 0; i < m_trialLocalOps.size(); ++i)
    discreteTrialLocalOps[i] = m_trialLocalOps[i].weakForm();

  // Calculate the inverse mass matrices
  typedef DiscreteSparseBoundaryOperator<ResultType> SparseOp;
  shared_ptr<const SparseOp> discreteTestId, discreteTrialId;
  shared_ptr<RealSparseMatrix> testInverse, trialInverse;
  if (!discreteTestLocalOps.empty()) {
    BoundaryOp testId = identityOperator(
        // We don't need a persistent shared_ptr since identityOperator
        // will go out of scope at the end of this function anyway.
        // All we need is a weak form.
        auxContext, m_integralOp.dualToRange(), m_integralOp.dualToRange(),
        m_integralOp.dualToRange(), "(" + this->label() + ")_test_id");
    discreteTestId = dynamic_pointer_cast<const SparseOp>(testId.weakForm());
    if (!discreteTestId)
      throw std::runtime_error(
          "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
          "identity operator must be represented by a sparse matrix");

    // NOTE: here we form the explicit inverse of the mass matrix. Of
    // course it is not the right way to do it and we'd be better off using
    // Cholesky decomposition and triangular solves. However, implementing
    // this for sparse matrices would be relatively hard. Moreover, in
    // practice the mass matrix will be block diagonal (since the internal
    // domain and dual to range spaces will be discontinuous), so we only
    // form explicit inverses of very small (e.g. 3 x 3) blocks, and mass
    // matrices are usually well conditioned. So explicit inversion
    // shouldn't cause any significant loss of accuracy.
    testInverse = sparseInverse(*discreteTestId->sparseMatrix());
  }

  if (!discreteTrialLocalOps.empty()) {
    if (m_integralOp.domain() == m_integralOp.dualToRange() && discreteTestId) {
      discreteTrialId = discreteTestId;
      //            discreteTrialInvId = discreteTestInvId;
      trialInverse = testInverse;
    } else {
      BoundaryOp trialId = identityOperator(
          auxContext, m_integralOp.domain(), m_integralOp.domain(),
          m_integralOp.domain(), "(" + this->label() + ")_trial_id");
      discreteTrialId =
          dynamic_pointer_cast<const SparseOp>(trialId.weakForm());
      if (!discreteTrialId)
        throw std::runtime_error(
            "SyntheticIntegralOperator::SyntheticIntegralOperator(): "
            "identity operator must be represented by a sparse matrix");
      trialInverse = sparseInverse(*discreteTrialId->sparseMatrix());
    }
  }

  // Coalesce sparse operators
  discreteTestLocalOps =
      coalesceTestOperators(discreteTestLocalOps, testInverse);
  if (discreteTrialLocalOps.empty()) {
    if (symmetricMode)
      discreteTrialLocalOps = transposeTestOperators(
          discreteTestLocalOps, m_syntheseSymmetry & HERMITIAN);
  } else
    discreteTrialLocalOps =
        coalesceTrialOperators(discreteTrialLocalOps, trialInverse);

  // Now join all the pieces together
  shared_ptr<DiscreteLinOp> result;
  if (discreteTestLocalOps.empty()) {
    result = mul<ResultType>(discreteIntegralOp, discreteTrialLocalOps[0]);
    for (size_t i = 1; i < discreteTestLocalOps.size(); ++i)
      result =
          sum<ResultType>(result, mul<ResultType>(discreteIntegralOp,
                                                  discreteTrialLocalOps[i]));
  } else if (discreteTrialLocalOps.empty()) {
    result = mul<ResultType>(discreteTestLocalOps[0], discreteIntegralOp);
    for (size_t i = 1; i < discreteTestLocalOps.size(); ++i)
      result = sum<ResultType>(
          result, mul<ResultType>(discreteTestLocalOps[i], discreteIntegralOp));
  } else {
    result = mul<ResultType>(
        mul<ResultType>(discreteTestLocalOps[0], discreteIntegralOp),
        discreteTrialLocalOps[0]);
    for (size_t i = 1; i < discreteTestLocalOps.size(); ++i)
      result = sum<ResultType>(
          result, mul<ResultType>(mul<ResultType>(discreteTestLocalOps[i],
                                                  discreteIntegralOp),
                                  discreteTrialLocalOps[i]));
  }

  tbb::tick_count end = tbb::tick_count::now();

  if (verbose)
    std::cout << "Assembly of the weak form of operator '" << this->label()
              << "' took " << (end - start).seconds() << " s" << std::endl;
  return result;
}

template <typename BasisFunctionType, typename ResultType>
void SyntheticIntegralOperator<BasisFunctionType, ResultType>::
    getContextsForInternalAndAuxiliaryOperators(
        const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
        shared_ptr<const Context<BasisFunctionType, ResultType>>
            &internalContext,
        shared_ptr<const Context<BasisFunctionType, ResultType>> &auxContext) {
  typedef Context<BasisFunctionType, ResultType> Ctx;
  ParameterList parameters(context->globalParameterList());
  parameters.put("options.hmat.hMatAssemblyMode",
                 std::string("GlobalAssembly"));
  internalContext.reset(new Ctx(parameters));
  parameters.put("options.global.verbosityLevel", static_cast<int>(-5));
  auxContext.reset(new Ctx(parameters));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    SyntheticIntegralOperator);

} // namespace Bempp
