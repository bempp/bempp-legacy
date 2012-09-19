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

#ifndef bempp_numerical_quadrature_strategy_hpp
#define bempp_numerical_quadrature_strategy_hpp

#include "../common/common.hpp"

#include "../fiber/numerical_quadrature_strategy.hpp"
#include "../grid/geometry_factory.hpp"

namespace Bempp
{

using Fiber::AccuracyOptions;
using Fiber::AccuracyOptionsEx;

/** \ingroup weak_form_assembly
 *  \brief Numerical quadrature strategy.
 *
 *  A quadrature strategy provides functions constructing local assemblers used
 *  to discretize boundary operators and user-defined functions. A particular
 *  quadrature strategy determines how the integrals involved in this
 *  discretization are evaluated.
 *
 *  The local assemblers constructed by this class use numerical quadrature to
 *  evaluate the necessary integrals. Singular integrals are transformed into
 *  regular ones as described in S. Sauter, Ch. Schwab, "Boundary Element
 *  Methods" (2010). Quadrature accuracy can be influenced by parameters given
 *  during the construction. */
template <typename BasisFunctionType, typename ResultType>
class NumericalQuadratureStrategy :
        public Fiber::NumericalQuadratureStrategy<
        BasisFunctionType, ResultType, GeometryFactory>
{
private:
    typedef Fiber::NumericalQuadratureStrategy<
    BasisFunctionType, ResultType, GeometryFactory> Base;
public:
    /** \brief Construct a numerical quadrature strategy with default accuracy
     *  settings.
     *
     *  Calling this constructor is equivalent to calling the other constructor
     *  with \p accuracyOptions equal to <tt>AccuracyOptionsEx()</tt>. */
    NumericalQuadratureStrategy();

    /** \brief Construct a numerical quadrature strategy with prescribed
     *  accuracy settings.
     *
     *  The quadrature order is determined differently for different types of
     *  integrals:
     *
     *    <ol>
     *    <li> <b>Double integrals of the form
     *          \f[ \int_{\Gamma} \int_{\Sigma} f(x, y) \,
     *              d\Gamma(x) \, d\Sigma(y), \f]
     *      where \f$\Gamma\f$ and \f$\Sigma\f$ are two disjoint elements and
     *      \f$f(x, y)\f$ is a function regular for \f$x \in \Gamma\f$ and \f$y
     *      \in \Sigma\f$.</b> Such integrals are approximated by
     *          \f[ \sum_{i=1}^m \sum_{j=1}^n w_i^m w_j^n \, f(x_i^m, y_j^n), \f]
     *      where \f$x_i^m\f$ and \f$y_j^n\f$ are appropriate quadrature points
     *      and \f$w_i^m\f$ and \f$w_j^n\f$ are the corresponding quadrature
     *      weights. By default, these are chosen so that the order of accuracy
     *      of the quadrature of each element is equal to the maximum degree of
     *      the polynomials belonging to the basis attached to that element. In
     *      other words, the quadrature rule is chosen so that a function
     *      \f$f(x, y)\f$ being a product of two polynomials, \f$u(x)\f$ and
     *      \f$v(y)\f$, with degrees equal to the orders of the bases attached
     *      to elements \f$\Gamma\f$ and \f$\Sigma\f$ would be integrated
     *      exactly. For instance, for a pair of elements endowed with linear
     *      bases, single-point quadrature is by default used on both elements.
     *
     *      This default integration order is often insufficient. It can be
     *      increased by calling one ofthe overloads of the
     *      AccuracyOptionsEx::setDoubleRegular() function on the object passed
     *      subsequently to this constructor via the \p accuracyOptions
     *      parameter. The quadrature rule can be made dependent on the
     *      distance between the two elements; see the documentation of
     *      AccuracyOptionsEx for details.
     *
     *    <li> <b>Double integrals of the same form as above, but on
     *      pairs of elements \f$\Gamma\f$ and \f$\Sigma\f$ sharing at least a
     *      single point and with the function \f$f(x, y)\f$ having a
     *      singularity at \f$x = y\f$.</b> Such integrals are evaluated by first
     *      applying an appropriate coordinate transformation to remove the
     *      singularity of the integrand, as described in the book of Sauter
     *      and Schwab cited before, and then approximating the new integral
     *      with a tensor-product Gaussian quadrature rule with order of
     *      accuracy in each dimension choosen by default as \f$\max(p, q) +
     *      5\f$, where \f$p\f$ and \f$q\f$ are the orders of the bases
     *      attached to elements \f$\Gamma\f$ and \f$\Sigma\f$. The order of
     *      accuracy can be customized by calling the
     *      AccuracyOptionsEx::setDoubleSingular() function.
     *
     *    <li> <b>Integrals over single elements
     *          \f[ \int_{\Gamma} f(x) \, d\Gamma(x) \f]
     *      of regular functions \f$f(x)\f$.</b> They are evaluated using a
     *      Gaussian quadrature rule with order of accuracy taken by default as
     *      twice the order of the basis attached to the element \f$\Gamma\f$.
     *      This order of accuracy can be customized by calling the
     *      AccuracyOptionsEx::setSingleRegular() function.
     *    </ul>
     *
     *  \note For backward compatibility reasons, an object of type
     *  AccuracyOptions can be passed into this constructor instead of an
     *  AccuracyOptionsEx object. The use of the AccuracyOptions class is
     *  deprecated.
     */
    explicit NumericalQuadratureStrategy(
            const AccuracyOptionsEx& accuracyOptions);
};

} // namespace Bempp

#endif
