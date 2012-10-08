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

#ifndef bempp_discrete_boundary_operator_hpp
#define bempp_discrete_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include "transposition_mode.hpp"
#include "boost/enable_shared_from_this.hpp"

#include "../common/armadillo_fwd.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_LinearOpDefaultBase_decl.hpp>
#endif

namespace Bempp
{

/** \ingroup discrete_boundary_operators
 *  \brief Discrete boundary operator.
 *
 *  This class represents a discretised boundary operator, i.e. a matrix
 *  \f$\mathsf{L}\f$ with entries
 *  \f[ \mathsf{L}_{mn} \equiv \int_S \phi^*_m(x) \,[L\, \psi_n](x)\,
 *      \mathrm{d}S(x), \f]
 *  where \f$L\f$ is a boundary operator, \f$\{\phi_m\}_{m=1}^{M}\f$ are the
 *  basis functions spanning a *test space* and defined on a surface \f$S\f$,
 *  whereas \f$\{\psi_n\}_{n=1}^{N}\f$ are the basis functions of a *trial
 *  space*. The way the matrix is stored can be arbitrary. Concrete subclasses
 *  of this class implement specific storage methods.
 *
 *  If BEM++ is compiled with Trilinos, this class is derived from
 *  <tt>Thyra::LinearOpDefaultBase<ValueType></tt> and hence inherits all
 *  non-private member function of this class; see <a
 *  href="http://trilinos.sandia.gov/packages/docs/dev/packages/thyra/doc/html/classThyra_1_1LinearOpDefaultBase.html">Trilinos
 *  documentation</a> for full documentation of these functions.
 *
 *  If Trilinos is not available during compilation, a simple fallback interface
 *  is provided.
 *
 *  \tparam ValueType
 *    Type used to represent entries of the operator matrix.
 *    It can take the following values: \c float, \c
 *    double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 */
template <typename ValueType>
class DiscreteBoundaryOperator : public boost::enable_shared_from_this<DiscreteBoundaryOperator<ValueType> >
#ifdef WITH_TRILINOS
        , public Thyra::LinearOpDefaultBase<ValueType>
#endif
{
public:
    /** \brief Destructor. */
    virtual ~DiscreteBoundaryOperator() {}

#ifdef WITH_TRILINOS
#ifndef DOXYGEN
    // import the apply() member function from base class
    // before we overload it
    using Thyra::LinearOpDefaultBase<ValueType>::apply;
#else // DOXYGEN
    /** \brief Apply the linear operator to a multivector.
     *
     *  Set the elements of the multivector <tt>y_inout</tt> to
     *
     *  <tt>y_inout := alpha * trans(L) * x_in + beta * y_inout</tt>,
     *
     *  where \c L is the linear operator represented by this object.
     *
     *  \param[in] trans
     *    Determines whether what is applied is the "bare" operator, its
     *    transpose, conjugate or conjugate transpose.
     *  \param[in] x_in
     *    The right-hand-side multivector.
     *  \param[in,out] y_inout
     *    The target multivector being transformed. When <tt>beta == 0.0</tt>,
     *    this multivector can have uninitialized elements.
     *  \param[in] alpha
     *    Scalar multiplying this operator.
     *  \param[in] beta
     *    The multiplier for the target multivector \p y_inout.
     *
     *  This overload of the apply() member function is only available if
     *  the library was compiled with Trilinos. */
    void apply(const Thyra::EOpTransp M_trans,
               const Thyra::MultiVectorBase<ValueType>& x_in,
               const Thyra::Ptr<Thyra::MultiVectorBase<ValueType> >& y_inout,
               const ValueType alpha,
               const ValueType beta) const;
#endif // DOXYGEN
#endif // WITH_TRILINOS

    /** \brief Apply the operator to a vector.
     *
     *  Set the elements of the column vector \p y_inout to
     *
     *  <tt>y_inout := alpha * trans(L) * x_in + beta * y_inout</tt>,
     *
     *  where \c L is the linear operator represented by this object.
     *
     *  \param[in] trans
     *    Determines whether what is applied is the "bare" operator, its
     *    transpose, conjugate or conjugate transpose.
     *  \param[in] x_in
     *    The right-hand-side multivector.
     *  \param[in,out] y_inout
     *    The target multivector being transformed. When <tt>beta == 0.0</tt>,
     *    this multivector can have uninitialized elements.
     *  \param[in] alpha
     *    Scalar multiplying this operator.
     *  \param[in] beta
     *    The multiplier for the target multivector \p y_inout.
     *
     *  This overload is always available, even if the library was compiled
     *  without Trilinos.
     */
    void apply(const TranspositionMode trans,
               const arma::Col<ValueType>& x_in,
               arma::Col<ValueType>& y_inout,
               const ValueType alpha,
               const ValueType beta) const {
        bool transposed = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);
        if (x_in.n_rows != (transposed ? rowCount() : columnCount()))
            throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                        "vector x_in has invalid length");
        if (y_inout.n_rows != (transposed ? columnCount() : rowCount()))
            throw std::invalid_argument("DiscreteBoundaryOperator::apply(): "
                                        "vector y_inout has invalid length");

        applyBuiltInImpl(trans, x_in, y_inout, alpha, beta);
    }

    /** \brief Return a representation that can be casted to a DiscreteAcaBoundaryOperator
      *
      * The conversion only succeeds if all members if the DiscreteOperator itself can be
      * cast to a DiscreteAcaBoundaryOperator or if it is a linear combination of
      * DiscreteOperators that can be cast to type DiscreteAcaBoundaryOperator.
      * Operator compositions are not yet supported by this function.
      *
      * \param[in] eps
      * Accuracy tolerance for H-Matrix addition.
      * \param[in] maximumRank
      * maximum rank of blocks to be considered low rank in the case of H-Matrix addition.
      * \returns A pointer to a DiscreteBoundaryOperator object, which is castable to
      * DiscreteAcaBoundaryOperator.
      */
    virtual shared_ptr<const DiscreteBoundaryOperator>
    asDiscreteAcaBoundaryOperator(double eps=1E-4, int maximumRank=50) const = 0;

    /** \brief Write a textual representation of the operator to standard output.
     *
     *  The default implementation prints the matrix returned by asMatrix().
     */
    virtual void dump() const;

    /** \brief Matrix representation of the operator.

    The default implementation is slow and should be overridden where possible. */
    virtual arma::Mat<ValueType> asMatrix() const;

    /** \brief Number of rows of the operator. */
    virtual unsigned int rowCount() const = 0;

    /** \brief Number of columns of the operator. */
    virtual unsigned int columnCount() const = 0;

    /** \brief Add a subblock of this operator to a matrix.
     *
     *  Perform the operation <tt>block += alpha * L[rows, cols]</tt>,
     *  where \p block is a matrix, \p alpha a scalar and <tt>L[rows, cols]</tt>
     *  a subblock of this discrete linear operator.
     *
     *  \param[in] rows  %Vector of row indices.
     *  \param[in] cols  %Vector of column indices.
     *  \param[in] alpha Multiplier.
     *  \param[in,out] block
     *    On entry, matrix of size (<tt>rows.size()</tt>, <tt>cols.size()</tt>).
     *    On exit, each element (\e i, \e j) of this matrix will be augmented by
     *    the element (<tt>rows[i]</tt>, <tt>cols[j]</tt>) of this operator,
     *    multiplied by \p alpha.
     *
     *  Row and column indices may be unsorted.
     *
     * This method need not be supported by all subclasses. It is mainly
     * intended for internal use during weak-form assembly. */
    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const = 0;

#ifdef WITH_TRILINOS
protected:
    virtual void applyImpl(
            const Thyra::EOpTransp M_trans,
            const Thyra::MultiVectorBase<ValueType> &X_in,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> > &Y_inout,
            const ValueType alpha,
            const ValueType beta) const;
#endif

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const = 0;
};

/** \brief Unary plus: return a copy of the argument. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator+(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \brief Return a shared pointer to a DiscreteBoundaryOperator representing
 *  the operand multiplied by -1.
 *
 *  An exception is thrown if \p op is null. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator-(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \brief Return a shared pointer to a DiscreteBoundaryOperator representing
 *  the sum of the operands.
 *
 *  An exception is thrown if either \p op1 or \p op2 is null. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator+(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2);

/** \brief Return a shared pointer to a DiscreteBoundaryOperator representing
 *  the difference of the operands.
 *
 *  An exception is thrown if either \p op1 or \p op2 is null. */
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator-(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op1,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op2);

/** \brief Return a shared pointer to a DiscreteBoundaryOperator representing
 *  the operator <tt>*op</tt> multiplied by \p scalar.
 *
 *  An exception is thrown if \p op is null. */
template <typename ValueType, typename ScalarType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator*(
        ScalarType scalar,
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op);

/** \brief Return a shared pointer to a DiscreteBoundaryOperator representing
 *  the operator <tt>*op</tt> multiplied by \p scalar.
 *
 *  An exception is thrown if \p op is null. */
template <typename ValueType, typename ScalarType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator*(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        ScalarType scalar);

/** \brief Return a shared pointer to a DiscreteBoundaryOperator representing
 *  the operator <tt>*op</tt> divided by \p scalar.
 *
 *  An exception is thrown if \p op is null. */
template <typename ValueType, typename ScalarType>
shared_ptr<const DiscreteBoundaryOperator<ValueType> > operator/(
        const shared_ptr<const DiscreteBoundaryOperator<ValueType> >& op,
        ScalarType scalar);

} // namespace Bempp

#endif
