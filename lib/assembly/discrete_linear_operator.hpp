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

#include "config_trilinos.hpp"

#ifndef bempp_discrete_linear_operator_hpp
#define bempp_discrete_linear_operator_hpp

#include "../common/common.hpp"

#include "transposition_mode.hpp"

#include "../common/armadillo_fwd.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_LinearOpDefaultBase_decl.hpp>
#endif

namespace Bempp
{

/** \ingroup assembly
 *  \brief Discrete linear operator.
 *
 *  This class represents a discretised linear operator, i.e. a matrix
 *  \f$\mathsf{L}\f$ with entries
 *  \f[ \mathsf{L}_{mn} \equiv \int_S \phi^*_m(x) \,[L\, \psi_n](x)\,
 *      \mathrm{d}S(x), \f]
 *  where \f$L\f$ is a linear operator acting on functions defined on a surface
 *  \f$S\f$, \f$\{\phi_m\}_{m=1}^{M}\f$ is a basis of the test function space and
 *  \f$\{\psi_n\}_{n=1}^{N}\f$ is a basis of the trial function space.
 *  The way the matrix is stored can be arbitrary. Concrete subclasses of
 *  this class implement specific storage methods.
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
class DiscreteLinearOperator
#ifdef WITH_TRILINOS
        : public Thyra::LinearOpDefaultBase<ValueType>
#endif
{
public:
    /** \brief Destructor. */
    virtual ~DiscreteLinearOperator() {}

#ifdef WITH_TRILINOS
    // import the apply() member function from base class
    // before we overload it
    using Thyra::LinearOpDefaultBase<ValueType>::apply;
#endif

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
        applyBuiltInImpl(trans, x_in, y_inout, alpha, beta);
    }

    /** \brief Write a textual representation of the operator to standard output. */
    virtual void dump() const = 0;

    /** \brief Matrix representation of the operator.

    This method need not be supported by all subclasses. */
    virtual arma::Mat<ValueType> asMatrix() const = 0;

    /** \brief Number of rows of the operator. */
    virtual unsigned int rowCount() const = 0;

    /** \brief Number of columns of the operator. */
    virtual unsigned int columnCount() const = 0;

    /** \brief Perform the operation <tt>block += alpha * L[rows, cols]</tt>,
      where \p alpha is a scalar and <tt>L[rows, cols]</tt> is a subblock of this
      discrete linear operator.

      \param[in] rows  %Vector of row indices.
      \param[in] cols  %Vector of col indices.
      \param[in] alpha scalar multiplier for the block
      \param[in,out] block
        On entry, matrix of size (<tt>rows.size()</tt>, <tt>cols.size()</tt>).
        On exit, each element (\e i, \e j) of this matrix will be augmented by
        the element (<tt>rows[i]</tt>, <tt>cols[j]</tt>) of this operator,
        multiplied by \p alpha.

      Row and column indices may be unsorted.

      This method need not be supported by all subclasses. It is mainly intended
      for internal use during weak-form assembly. */
    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const = 0;

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const = 0;
};

} // namespace Bempp

#endif
