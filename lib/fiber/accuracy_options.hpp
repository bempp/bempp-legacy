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

#ifndef fiber_accuracy_options_hpp
#define fiber_accuracy_options_hpp

#include "../common/common.hpp"

#include "quadrature_options.hpp"

#include <limits>
#include <utility>
#include <vector>

namespace Fiber
{

/** \brief Old-style options controlling quadrature accuracy.
 *
 *  \deprecated Use AccuracyOptionsEx instead. */
struct AccuracyOptions
{
public:
    /** \brief Options controlling integration of regular functions
     *  on single elements. */
    QuadratureOptions singleRegular;
    /** \brief Options controlling integration of regular functions
     *  on pairs of elements. */
    QuadratureOptions doubleRegular;
    /** \brief Options controlling integration of singular functions
     *  on pairs of elements. */
    QuadratureOptions doubleSingular;
};

/** \brief New-style options controlling quadrature accuracy. */
class AccuracyOptionsEx
{
public:
    /** \brief Constructor.
     *
     *  Create an AccuracyOptionsEx object representing default quadrature
     *  accuracy settings. */
    AccuracyOptionsEx();

    /** \brief Constructor.
     *
     *  Convert an old-style AccuracyOptions object to an equivalent new-style
     *  AccuracyOptionsEx object.
     *
     *  \note This constructor is non-explicit on purpose (to allow an implicit
     *  conversion from an old-style to a new-style options object). */
    AccuracyOptionsEx(const AccuracyOptions& oldStyleOpts);

    /** \brief Return the options controlling integration of regular functions
     *  on single elements. */
    const QuadratureOptions& singleRegular() const;

    /** \brief Return the options controlling integration, on single elements,
     *  of functions that have singularities outside these elements.
     *
     *  This function is called, in particular, to determine the quadrature order
     *  used in the approximation of integrals of the form
     *
     *  \f[ \int_E G(x, y) f(y) \, \mathrm{d}(y), \f]
     *
     *  where \f$E\f$ is an element, \f$f(y)\f$ a regular function and \f$G(x,
     *  y)\f$ a function singular at \f$x = y\f$. The \p normalizedDistance
     *  argument should be set to an estimate of the minimum distance between
     *  the point \f$x\f$ and the element \f$E\f$, divided by the size of
     *  \f$E\f$ i.e. \f$\min_{y \in E} \lvert x - y \rvert / \lvert E\rvert\f$.
     */
    const QuadratureOptions& singleRegular(double normalizedDistance) const;

    /** \brief Set the options controlling integration of functions
     *  on single elements.
     *
     *  Use this overload of setSingleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  functions on single elements (or its increase \f$\Delta o\f$ over the
     *  default level if \p relativeToDefault is set to \c true) should be set
     *  to \p accuracyOrder, regardless of whether the integrated function has
     *  a singularity near the element or not. */
    void setSingleRegular(int accuracyOrder, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of functions
     *  on single elements.
     *
     *  Use this overload of setSingleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  functions on single elements (or its increase \f$\Delta o\f$ over the
     *  default level if \p relativeToDefault is set to \c true) should be
     *  chosen as follows, depending on the estimated normalized distance \em d
     *  of the element to the nearest singularity of the function being
     *  integrated:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2.
     *
     *  The normalized distance \em d is defined in the documentation of the
     *  singleRegular() method.
     */
    void setSingleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          int accuracyOrder2, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of functions
     *  on single elements.
     *
     *  Use this overload of setSingleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  functions on single elements (or its increase \f$\Delta o\f$ over the
     *  default level if \p relativeToDefault is set to \c true) should be
     *  chosen as follows, depending on the estimated normalized distance \em d
     *  of the element to the nearest singularity of the function being
     *  integrated:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - if \em d <= \p maxNormalizedDistance2:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder3.
     *
     *  The normalized distance \em d is defined in the documentation of the
     *  singleRegular() method.
     */
    void setSingleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          double maxNormalizedDistance2, int accuracyOrder2,
                          int accuracyOrder3, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of functions
     *  on single elements.
     *
     *  Use this overload of setSingleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  functions on single elements (or its increase \f$\Delta o\f$ over the
     *  default level if \p relativeToDefault is set to \c true) should be
     *  chosen as follows, depending on the estimated normalized distance \em d
     *  of the element to the nearest singularity of the function being
     *  integrated:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - if \em d <= \p maxNormalizedDistance2:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2
     *  - if \em d <= \p maxNormalizedDistance3:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder3
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder4.
     *
     *  The normalized distance \em d is defined in the documentation of the
     *  singleRegular() method.
     */
    void setSingleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          double maxNormalizedDistance2, int accuracyOrder2,
                          double maxNormalizedDistance3, int accuracyOrder3,
                          int accuracyOrder4, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of functions
     *  on single elements.
     *
     *  Use this overload of setSingleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  functions on single elements (or its increase \f$\Delta o\f$ over the
     *  default level if \p relativeToDefault is set to \c true) should be
     *  chosen as follows, depending on the estimated normalized distance \em d
     *  of the element to the nearest singularity of the function being
     *  integrated:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - if \em d <= \p maxNormalizedDistance2:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2
     *  - if \em d <= \p maxNormalizedDistance3:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder3
     *  - if \em d <= \p maxNormalizedDistance4:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder4
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder5.
     *
     *  The normalized distance \em d is defined in the documentation of the
     *  singleRegular() method.
     */
    void setSingleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          double maxNormalizedDistance2, int accuracyOrder2,
                          double maxNormalizedDistance3, int accuracyOrder3,
                          double maxNormalizedDistance4, int accuracyOrder4,
                          int accuracyOrder5, bool relativeToDefault = true);

    void setSingleRegular(const std::vector<double>& maxNormalizedDistances,
                          const std::vector<int>& accuracyOrders,
                          bool relativeToDefault = true);

    /** \brief Return the options controlling integration of regular functions
     *  on pairs of elements.
     *
     *  The quadrature rule used to approximate this class of integrals may
     *  depend on the distance between the two the elements. To retrieve the
     *  quadrature rule that should be used to integrate a function defined on
     *  a given pair of elements, set \p normalizedDistance to the
     *  distance between the centres of the elements divided by the size of the
     *  larger of the two elements. */
    const QuadratureOptions& doubleRegular(double normalizedDistance) const;

    /** \brief Set the options controlling integration of regular functions
     *  on pairs of elements.
     *
     *  Use this overload of setDoubleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  regular functions on pairs of elements (or its increase \f$\Delta o\f$
     *  over the default level if \p relativeToDefault is set to \c true)
     *  should be set to \p accuracyOrder, regardless of the distance between
     *  the two elements. */
    void setDoubleRegular(int accuracyOrder, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of regular functions
     *  on pairs of elements.
     *
     *  Use this overload of setDoubleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  regular functions on pairs of elements (or its increase \f$\Delta o\f$
     *  over the default level if \p relativeToDefault is set to \c true)
     *  should be chosen as follows, depending on the normalized distance \em d
     *  of the two elements:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2.
     *
     *  The normalized distance \em d is the distance between the centres of
     *  the two elements divided by the size of the larger of the elements. */
    void setDoubleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          int accuracyOrder2, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of regular functions
     *  on pairs of elements.
     *
     *  Use this overload of setDoubleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  regular functions on pairs of elements (or its increase \f$\Delta o\f$
     *  over the default level if \p relativeToDefault is set to \c true)
     *  should be chosen as follows, depending on the normalized distance \em d
     *  of the two elements:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - if \em d <= \p maxNormalizedDistance2:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder3.
     *
     *  The normalized distance \em d is the distance between the centres of
     *  the two elements divided by the size of the larger of the elements. */
    void setDoubleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          double maxNormalizedDistance2, int accuracyOrder2,
                          int accuracyOrder3, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of regular functions
     *  on pairs of elements.
     *
     *  Use this overload of setDoubleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  regular functions on pairs of elements (or its increase \f$\Delta o\f$
     *  over the default level if \p relativeToDefault is set to \c true)
     *  should be chosen as follows, depending on the normalized distance \em d
     *  of the two elements:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - if \em d <= \p maxNormalizedDistance2:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2
     *  - if \em d <= \p maxNormalizedDistance3:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder3
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder4.
     *
     *  The normalized distance \em d is the distance between the centres of
     *  the two elements divided by the size of the larger of the elements. */
    void setDoubleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          double maxNormalizedDistance2, int accuracyOrder2,
                          double maxNormalizedDistance3, int accuracyOrder3,
                          int accuracyOrder4, bool relativeToDefault = true);

    /** \brief Set the options controlling integration of regular functions
     *  on pairs of elements.
     *
     *  Use this overload of setDoubleRegular() to specify that the order of
     *  accuracy \em o of the quadrature rule used to approximate integrals of
     *  regular functions on pairs of elements (or its increase \f$\Delta o\f$
     *  over the default level if \p relativeToDefault is set to \c true)
     *  should be chosen as follows, depending on the normalized distance \em d
     *  of the two elements:
     *
     *  - if \em d <= \p maxNormalizedDistance1:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder1
     *  - if \em d <= \p maxNormalizedDistance2:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder2
     *  - if \em d <= \p maxNormalizedDistance3:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder3
     *  - if \em d <= \p maxNormalizedDistance4:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder4
     *  - else:
     *    - set \em o or \f$\Delta o\f$ to accuracyOrder5.
     *
     *  The normalized distance \em d is the distance between the centres of
     *  the two elements divided by the size of the larger of the elements. */
    void setDoubleRegular(double maxNormalizedDistance1, int accuracyOrder1,
                          double maxNormalizedDistance2, int accuracyOrder2,
                          double maxNormalizedDistance3, int accuracyOrder3,
                          double maxNormalizedDistance4, int accuracyOrder4,
                          int accuracyOrder5, bool relativeToDefault = true);

    void setDoubleRegular(const std::vector<double>& maxNormalizedDistances,
                          const std::vector<int>& accuracyOrders,
                          bool relativeToDefault = true);

    /** \brief Return the options controlling integration of singular functions
     *  on pairs of elements. */
    const QuadratureOptions& doubleSingular() const;

    /** \brief Set the options controlling integration of singular functions
     *  on pairs of elements.
     *
     *  If \p relativeToDefault is set to \c false, \p accuracyOrder denotes
     *  the desired order of accuracy of the quadrature rule used to approximate
     *  integrals of singular functions on pairs of elements. Otherwise,
     *  \p accuracyOrder denotes the desired *increase* of the order of accuracy
     *  above the default level. */
    void setDoubleSingular(int accuracyOrder, bool relativeToDefault = true);

private:
    /** \cond PRIVATE */
    std::vector<std::pair<double, QuadratureOptions> > m_singleRegular;
    std::vector<std::pair<double, QuadratureOptions> > m_doubleRegular;
    QuadratureOptions m_doubleSingular;
    /** \endcond */
};

} // namespace Fiber

#endif
