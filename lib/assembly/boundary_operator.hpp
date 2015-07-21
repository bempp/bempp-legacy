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

#ifndef bempp_boundary_operator_hpp
#define bempp_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"
#include "transposition_mode.hpp"

#include <boost/mpl/set.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/weak_ptr.hpp>
#include <string>
#include <iostream>

#include <complex>

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Space;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType>
class AbstractBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class Context;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
/** \endcond */

/** \ingroup boundary_operators
 *  \brief Operator acting on functions defined on a surface.
 *
 *  A BoundaryOperator is a lightweight wrapper of a pair of shared pointers to
 *  an AbstractBoundaryOperator and a DiscreteBoundaryOperator representing the
 *  weak form of the former. The weak form is evaluated lazily, on the first
 *  call to weakForm(). The Context object passed to the constructor of the
 *  BoundaryOperator or to the initialize() function determines how this weak
 *  form is calculated.
 *
 *  \note Different threads should not share BoundaryOperator objects, since
 *  notably the weakForm() function is not thread-safe. Instead, each thread
 *  should hold its own copy of a BoundaryOperator (note that copying
 *  BoundaryOperators is cheap -- the copy constructor is shallow).
 *
 *  See the documentation of AbstractBoundaryOperator for the decription of the
 *  template parameters \p BasisFunctionType and \p ResultType. */
template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator {
public:
  /** \brief Construct an uninitialized BoundaryOperator. */
  BoundaryOperator();

  /** \brief Destructor. */
  virtual ~BoundaryOperator();

  /** \brief Construct and initialize a BoundaryOperator.
   *
   *  Equivalent to calling the initialize() function on a BoundaryOperator
   *  object created with the default constructor. See the documentation of
   *  initialize() for a description of the constructor's parameters.
   *
   *  \note User code should not need to invoke this constructor directly; it
   *  is more convenient to use non-member constructors supplied with
   *  particular AbstractBoundaryOperator subclasses
   *  (e.g. laplace3dSingleLayerBoundaryOperator(), identityOperator(),
   *  pseudoinverse() etc.), which construct an AbstractBoundaryOperator and
   *  wrap it in a BoundaryOperator in a single step.
   */
  BoundaryOperator(
      const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
      const shared_ptr<const AbstractBoundaryOperator<BasisFunctionType,
                                                      ResultType>> &abstractOp);

  /** \brief Initialize or reinitialize a BoundaryOperator.
   *
   *  \param[in] context
   *    Shared pointer to a Context object that will be used
   *    to build the weak form of \p abstractOp when necessary.
   *  \param[in] abstractOp
   *    Shared pointer to an AbstractBoundaryOperator that will be
   *    encapsulated in this BoundaryOperator.
   *
   *  An exception is thrown if either of these pointers is NULL.
   *
   *  The provided shared pointers are stored in internal variables. In
   *  addition, any stored pointer to the weak form of the abstract boundary
   *  operator is invalidated. */
  void initialize(
      const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
      const shared_ptr<const AbstractBoundaryOperator<BasisFunctionType,
                                                      ResultType>> &abstractOp);

  /** \brief Uninitialize the BoundaryOperator.
   *
   *  This function resets the internal shared pointers to the abstract
   *  boundary operator and its weak form to NULL. */
  void uninitialize();

  /** \brief Return true if the BoundaryOperator has been initialized, false
   *  otherwise. */
  bool isInitialized() const;

  /** \brief Return a shared pointer to the encapsulated abstract boundary
   *  operator. */
  shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType>>
  abstractOperator() const;

  /** \brief Return a shared pointer to the stored Context object. */
  shared_ptr<const Context<BasisFunctionType, ResultType>> context() const;

  /** \brief Return a shared pointer to the weak form of the encapsulated
   *  abstract boundary operator.
   *
   *  An exception is thrown if this function is called on an uninitialized
   *  BoundaryOperator. If update is set to true then the weak form
   *  is recomputed.*/
  shared_ptr<const DiscreteBoundaryOperator<ResultType>> weakForm(bool update=false) const;

  /** \brief Return a shared pointer to the domain of the encapsulated
   *  abstract boundary operator.
   *
   *  A null pointer is returned if this function is called on an
   *  uninitialized BoundaryOperator. */
  shared_ptr<const Space<BasisFunctionType>> domain() const;

  /** \brief Return a shared pointer to the range of the encapsulated abstract
   *  boundary operator.
   *
   *  A null pointer is returned if this function is called on an
   *  uninitialized BoundaryOperator. */
  shared_ptr<const Space<BasisFunctionType>> range() const;

  /** \brief Return a shared pointer to the space dual to the range of the
   *  encapsulated abstract boundary operator.
   *
   *  A null pointer is returned if this function is called on an
   *  uninitialized BoundaryOperator. */
  shared_ptr<const Space<BasisFunctionType>> dualToRange() const;

  /** \brief Return the label of this BoundaryOperator. */
  std::string label() const;

  /** \brief Return true if the BoundaryOperator should prevent its weak form
   *  from being destroyed.
   *
   *  If isWeakFormHeld() returns true, the BoundaryOperator stores a
   *  "strong" shared pointer to its discrete weak form once it is assembled
   *  for the first time, so that the weak form is held in memory at least as
   *  long as the BoundaryOperator object is in scope and its uninitialize()
   *  member function is not called.
   *
   *  Otherwise only a weak pointer to the weak form is stored; thus, the
   *  weak form remains in memory only as long as a shared pointer to it is
   *  held somewhere else in the program.
   *
   *  By default, isWeakFormHeld() returns true; call holdWeakForm(false) to
   *  change it. */
  bool isWeakFormHeld() const;

  /** \brief Specify whether the BoundaryOperator should prevent its weak form
   *  from being destroyed.
   *
   *  See isWeakFormHeld() for more information. */
  void holdWeakForm(bool value);

  /** \brief Act on a GridFunction.
   *
   *  This function sets \p y_inout to <tt>alpha * A * x_in + beta *
   *  y_inout</tt>, where \c A is the operator represented by this object.
   *
   *  The space of \p x_in must be identical with the domain of the
   *  encapsulated abstract boundary operator, whereas the space of \p y_inout
   *  and its dual must be identical with the range of the encapsulated
   *  abstract boundary operator and its dual; otherwise an exception is
   *  thrown. An exception is also thrown if the BoundaryOperator is
   *  uninitialized. */
  void apply(const TranspositionMode trans,
             const GridFunction<BasisFunctionType, ResultType> &x_in,
             GridFunction<BasisFunctionType, ResultType> &y_inout,
             ResultType alpha, ResultType beta) const;

private:
  /** \cond PRIVATE */
  shared_ptr<const Context<BasisFunctionType, ResultType>> m_context;
  shared_ptr<const AbstractBoundaryOperator<BasisFunctionType, ResultType>>
      m_abstractOp;
  bool m_holdWeakForm;
  typedef shared_ptr<const DiscreteBoundaryOperator<ResultType>>
      ConstWeakFormContainer;
  mutable shared_ptr<ConstWeakFormContainer> m_weakFormContainer;
  typedef boost::weak_ptr<const DiscreteBoundaryOperator<ResultType>>
      WeakConstWeakFormContainer;
  mutable shared_ptr<WeakConstWeakFormContainer> m_weakWeakFormContainer;
  /** \endcond */
};

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator identical to the operand. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator+(const BoundaryOperator<BasisFunctionType, ResultType> &op);

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator equal to the operand multiplied by -1.
 *
 * An exception is thrown if either of the operands is uninitialized. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator-(const BoundaryOperator<BasisFunctionType, ResultType> &op);

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator representing the sum of the operands.
 *
 * An exception is thrown if either of the operands is uninitialized. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator+(const BoundaryOperator<BasisFunctionType, ResultType> &op1,
          const BoundaryOperator<BasisFunctionType, ResultType> &op2);

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator representing the difference of the
 *operands.
 *
 *  if either of the operands is uninitialized. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator-(const BoundaryOperator<BasisFunctionType, ResultType> &op1,
          const BoundaryOperator<BasisFunctionType, ResultType> &op2);

// This type machinery is needed to disambiguate between this operator and
// the one taking a BoundaryOperator and a GridFunction
/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator representing the operator \p op multiplied
 *  by \p scalar.
 *
 * An exception is thrown if \p op is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
typename boost::enable_if<
    typename boost::mpl::has_key<
        boost::mpl::set<float, double, std::complex<float>,
                        std::complex<double>>,
        ScalarType>,
    BoundaryOperator<BasisFunctionType, ResultType>>::type
operator*(const BoundaryOperator<BasisFunctionType, ResultType> &op,
          const ScalarType &scalar);

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator representing the operator \p op multiplied
 *  by \p scalar.
 *
 * An exception is thrown if \p op is uninitialized. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType>
operator*(const ScalarType &scalar,
          const BoundaryOperator<BasisFunctionType, ResultType> &op);

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator representing the operator \p op divided
 *  by \p scalar.
 *
 * \note \p scalar must not be zero. */
template <typename BasisFunctionType, typename ResultType, typename ScalarType>
BoundaryOperator<BasisFunctionType, ResultType>
operator/(const BoundaryOperator<BasisFunctionType, ResultType> &op,
          const ScalarType &scalar);

/** \relates BoundaryOperator
 *  \brief Act with a BoundaryOperator on a GridFunction.
 *
 *  This function returns the GridFunction obtained by acting with the operator
 *  \p op on the grid function \p fun. It is equivalent to calling
 *  \code
 *  op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
 *  return result;
 *  \endcode
 *  on GridFunction \p result with space and dual space compatible with
 *  the range and dual to range of \p op. */
template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
operator*(const BoundaryOperator<BasisFunctionType, ResultType> &op,
          const GridFunction<BasisFunctionType, ResultType> &fun);

/** \relates BoundaryOperator
 *  \brief Return a BoundaryOperator representing the product of the operands
 *  (<tt>op1 * op2</tt>).
 *
 *  An exception is thrown if any of the operands is uninitialized. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
operator*(const BoundaryOperator<BasisFunctionType, ResultType> &op1,
          const BoundaryOperator<BasisFunctionType, ResultType> &op2);

/** \relates BoundaryOperator
 *  \brief Return the adjoint of a BoundaryOperator.
 *
 *  An exception is thrown if the operand is uninitialized. This
 *  method tries to figure out the range space using a simple
 *  heuristic.
 *
 *  \see AdjointAbstractBoundaryOperator for the definition of the
 *  adjoint operator in BEM++. Note that this operator is only defined
 *  for real BasisFunctionType.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
adjoint(const BoundaryOperator<BasisFunctionType, ResultType> &op);

/** \relates BoundaryOperator
 *  \brief Return the adjoint of a BoundaryOperator.
 *
 *  An exception is thrown if the operand is uninitialized. The
 *  range space needs to be given explicitly.
 *
 *  \see AdjointAbstractBoundaryOperator for the definition of the
 *  adjoint operator in BEM++. Note that this operator is only defined
 *  for real BasisFunctionType.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
adjoint(const BoundaryOperator<BasisFunctionType, ResultType> &op,
        const shared_ptr<const Space<BasisFunctionType>> &range);

/** \relates BoundaryOperator
 *  \brief Check whether a BoundaryOperator object is initialized.
 *
 *  This function checks whether the BoundaryOperator \p op is initialized. If
 *  so, it returns a reference to \p op, otherwise it throws a
 *  <tt>std::invalid_argument</tt> exception with message \p message (or a
 *  default message if \p message is left empty). */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> &
throwIfUninitialized(BoundaryOperator<BasisFunctionType, ResultType> &op,
                     std::string message = "");

/** \relates BoundaryOperator
 *
 *  \overload */
template <typename BasisFunctionType, typename ResultType>
const BoundaryOperator<BasisFunctionType, ResultType> &
throwIfUninitialized(const BoundaryOperator<BasisFunctionType, ResultType> &op,
                     std::string message = "");

} // namespace Bempp

#endif
