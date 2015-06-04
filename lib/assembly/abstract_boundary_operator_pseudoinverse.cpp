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

#include "abstract_boundary_operator_pseudoinverse.hpp"

#include "context.hpp"
#include "discrete_inverse_sparse_boundary_operator.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "discrete_sparse_boundary_operator.hpp"

#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "../common/eigen_support.hpp"

#include "../assembly/discrete_inverse_sparse_boundary_operator.hpp"
#include "../assembly/discrete_sparse_boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator_composition.hpp"

#include <boost/make_shared.hpp>
#include <tbb/tick_count.h>

namespace Bempp {

////////////////////////////////////////////////////////////////////////////////
// AbstractBoundaryOperatorPseudoinverse

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>::
    AbstractBoundaryOperatorPseudoinverse(
        // TODO: add a solver argument specifying how to calculate the
        // pseudoinv.
        const BoundaryOperator<BasisFunctionType, ResultType> &operatorToInvert)
    : Base(operatorToInvert.range(), operatorToInvert.domain(),
           operatorToInvert.domain() == operatorToInvert.range()
               ? operatorToInvert.dualToRange()
               :
               // assume that domain == dualToRange, we'll verify it
               // in the body of the constructor
               operatorToInvert.range(),
           "pinv(" + operatorToInvert.label() + ")",
           throwIfUninitialized(operatorToInvert,
                                "AbstractBoundaryOperatorPseudoinverse::"
                                "AbstractBoundaryOperatorPseudoinverse(): "
                                "the boundary operator to be inverted must be "
                                "initialized")
               .abstractOperator()
               ->symmetry()),
      m_operator(operatorToInvert),
      m_id(boost::make_shared<AbstractBoundaryOperatorPseudoinverseId<
          BasisFunctionType, ResultType>>(operatorToInvert)) {
  if (operatorToInvert.domain() != operatorToInvert.range() &&
      operatorToInvert.domain() != operatorToInvert.dualToRange())
    throw std::runtime_error(
        "AbstractBoundaryOperatorPseudoinverse::"
        "AbstractBoundaryOperatorPseudoinverse(): "
        "Dual to the domain of the operator to invert cannot be determined "
        "since the domain is different from both "
        "the range and the space dual to its range");
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>::
    AbstractBoundaryOperatorPseudoinverse(
        // TODO: add a solver argument specifying how to calculate the
        // pseudoinv.
        const BoundaryOperator<BasisFunctionType, ResultType> &operatorToInvert,
        const shared_ptr<const Space<BasisFunctionType>> &dualToRange)
    : Base(operatorToInvert.range(), operatorToInvert.domain(), dualToRange,
           "pinv(" + operatorToInvert.label() + ")",
           throwIfUninitialized(operatorToInvert,
                                "AbstractBoundaryOperatorPseudoinverse::"
                                "AbstractBoundaryOperatorPseudoinverse(): "
                                "the boundary operator to be inverted must be "
                                "initialized")
               .abstractOperator()
               ->symmetry()),
      m_operator(operatorToInvert),
      m_id(boost::make_shared<AbstractBoundaryOperatorPseudoinverseId<
          BasisFunctionType, ResultType>>(operatorToInvert)) {}

template <typename BasisFunctionType, typename ResultType>
bool AbstractBoundaryOperatorPseudoinverse<BasisFunctionType,
                                           ResultType>::isLocal() const {
  return m_operator.abstractOperator()->isLocal();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>::id()
    const {
  return m_id;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>::
    assembleWeakFormImpl(
        const Context<BasisFunctionType, ResultType> &context) const {
  bool verbose =
      (context.assemblyOptions().verbosityLevel() >= VerbosityLevel::DEFAULT);
  shared_ptr<const DiscreteBoundaryOperator<ResultType>> wrappedDiscreteOp =
      m_operator.weakForm();

  if (verbose)
    std::cout << "Calculating the (pseudo)inverse of operator '"
              << m_operator.label() << "'..." << std::endl;

  tbb::tick_count start = tbb::tick_count::now();
  shared_ptr<DiscreteBoundaryOperator<ResultType>> result;
  if (shared_ptr<const DiscreteSparseBoundaryOperator<ResultType>>
          wrappedSparseOp = boost::dynamic_pointer_cast<
              const DiscreteSparseBoundaryOperator<ResultType>>(
              wrappedDiscreteOp))
    result = assembleWeakFormForSparseOperator(context, wrappedSparseOp);
  else if (shared_ptr<const DiscreteDenseBoundaryOperator<ResultType>>
               wrappedDenseOp = boost::dynamic_pointer_cast<
                   const DiscreteDenseBoundaryOperator<ResultType>>(
                   wrappedDiscreteOp))
    result = assembleWeakFormForDenseOperator(context, wrappedDenseOp);
  else
    throw std::runtime_error(
        "AbstractBoundaryOperatorPseudoinverse::assembleWeakFormImpl(): "
        "Currently only elementary boundary operators stored as sparse "
        "or dense matrices can be inverted");
  tbb::tick_count end = tbb::tick_count::now();

  if (verbose)
    std::cout << "Calculation of the (pseudo)inverse of operator '"
              << m_operator.label() << "' took " << (end - start).seconds()
              << " s" << std::endl;
  return result;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>::
    assembleWeakFormForSparseOperator(
        const Context<BasisFunctionType, ResultType> &context,
        const shared_ptr<const DiscreteSparseBoundaryOperator<ResultType>>
            &wrappedDiscreteOp) const {
  typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
  typedef DiscreteSparseBoundaryOperator<ResultType> DiscreteSparseOp;
  typedef DiscreteInverseSparseBoundaryOperator<ResultType>
      DiscreteInverseSparseOp;

  shared_ptr<const RealSparseMatrix> matrix = wrappedDiscreteOp->sparseMatrix();

  const int rowCount = matrix->rows();
  const int colCount = matrix->cols();
  if (rowCount == colCount) {
    // Square matrix; construct M^{-1}
    return boost::make_shared<DiscreteInverseSparseOp>(
        matrix, m_operator.abstractOperator()->symmetry());
  } else {
    // Construct the discrete operator representing M^H
    shared_ptr<DiscreteOp> transposeOp = boost::make_shared<DiscreteSparseOp>(
        matrix, NO_SYMMETRY, CONJUGATE_TRANSPOSE);

    // Construct the discrete operator representing the smaller of
    // (M^H M)^{-1} and (M M^H)^{-1}
    shared_ptr<RealSparseMatrix> productMatrix;

    if (rowCount > colCount) {
      // Tall matrix (overdetermined least-square problem);
      // construct (M^H M)^{-1} M^H
      productMatrix.reset(
          new RealSparseMatrix(matrix->transpose().eval() * (*matrix)));

      shared_ptr<DiscreteOp> productInverseOp =
          boost::make_shared<DiscreteInverseSparseOp>(productMatrix, HERMITIAN);

      return boost::make_shared<
          DiscreteBoundaryOperatorComposition<ResultType>>(productInverseOp,
                                                           transposeOp);
    } else {
      // Wide matrix (underdetermined least-square problem);
      // construct M^H (M M^H)^{-1}

      productMatrix.reset(
          new RealSparseMatrix((*matrix) * (matrix->transpose().eval())));
      shared_ptr<DiscreteOp> productInverseOp =
          boost::make_shared<DiscreteInverseSparseOp>(productMatrix, HERMITIAN);
      return boost::make_shared<
          DiscreteBoundaryOperatorComposition<ResultType>>(transposeOp,
                                                           productInverseOp);
    }
  }
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType>>
AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>::
    assembleWeakFormForDenseOperator(
        const Context<BasisFunctionType, ResultType> &context,
        const shared_ptr<const DiscreteDenseBoundaryOperator<ResultType>>
            &wrappedDiscreteOp) const {
  typedef DiscreteDenseBoundaryOperator<ResultType> DiscreteDenseLinOp;

  if (wrappedDiscreteOp->rowCount() == wrappedDiscreteOp->columnCount())
    // TODO: store an LU decomposition instead of the inverse matrix.
    return shared_ptr<DiscreteBoundaryOperator<ResultType>>(
        new DiscreteDenseLinOp(
            (wrappedDiscreteOp->asMatrix()).inverse().eval()));
  else
    // compute and store pseudoinverse
    return boost::make_shared<DiscreteDenseLinOp>(
        eigenMatPinv(wrappedDiscreteOp->asMatrix()));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> pseudoinverse(
    const BoundaryOperator<BasisFunctionType, ResultType> &boundaryOp) {
  typedef AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>
      Pinv;
  return BoundaryOperator<BasisFunctionType, ResultType>(
      boundaryOp.context(), boost::make_shared<Pinv>(boundaryOp));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
pseudoinverse(const BoundaryOperator<BasisFunctionType, ResultType> &boundaryOp,
              const shared_ptr<const Space<BasisFunctionType>> &dualToRange) {
  typedef AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>
      Pinv;
  return BoundaryOperator<BasisFunctionType, ResultType>(
      boundaryOp.context(), boost::make_shared<Pinv>(boundaryOp, dualToRange));
}

BEMPP_GCC_DIAG_OFF(deprecated - declarations);

////////////////////////////////////////////////////////////////////////////////
// AbstractBoundaryOperatorPseudoinverseId

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorPseudoinverseId<BasisFunctionType, ResultType>::
    AbstractBoundaryOperatorPseudoinverseId(
        const BoundaryOperator<BasisFunctionType, ResultType> &operatorToInvert)
    : m_operatorToInvertId(operatorToInvert.abstractOperator()->id()) {}

template <typename BasisFunctionType, typename ResultType>
size_t AbstractBoundaryOperatorPseudoinverseId<BasisFunctionType,
                                               ResultType>::hash() const {
  typedef AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>
      OperatorType;
  size_t result = tbb::tbb_hasher(typeid(OperatorType).name());
  tbb_hash_combine(result, m_operatorToInvertId->hash());
  return result;
}

template <typename BasisFunctionType, typename ResultType>
void AbstractBoundaryOperatorPseudoinverseId<BasisFunctionType,
                                             ResultType>::dump() const {
  typedef AbstractBoundaryOperatorPseudoinverse<BasisFunctionType, ResultType>
      OperatorType;
  std::cout << typeid(OperatorType).name() << " ";
  m_operatorToInvertId->dump();
}

template <typename BasisFunctionType, typename ResultType>
bool AbstractBoundaryOperatorPseudoinverseId<BasisFunctionType, ResultType>::
    isEqual(const AbstractBoundaryOperatorId &other) const {
  // dynamic_cast won't suffice since we want to make sure both objects
  // are of exactly the same type (dynamic_cast would succeed for a subclass)
  if (typeid(other) == typeid(*this)) {
    const AbstractBoundaryOperatorPseudoinverseId &otherCompatible =
        static_cast<const AbstractBoundaryOperatorPseudoinverseId &>(other);
    return (*m_operatorToInvertId == *(otherCompatible.m_operatorToInvertId));
  } else
    return false;
}

BEMPP_GCC_DIAG_ON(deprecated - declarations);

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT)                       \
  template BoundaryOperator<BASIS, RESULT> pseudoinverse(                      \
      const BoundaryOperator<BASIS, RESULT> &);                                \
  template BoundaryOperator<BASIS, RESULT> pseudoinverse(                      \
      const BoundaryOperator<BASIS, RESULT> &,                                 \
      const shared_ptr<const Space<BASIS>> &);
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    AbstractBoundaryOperatorPseudoinverse);

} // namespace Bempp
