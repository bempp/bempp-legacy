// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_function_hpp
#define fiber_function_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"

#include "../common/armadillo_fwd.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename CoordinateType> class GeometricalData;
/** \endcond */

/** \brief %Function intended to be evaluated on a boundary-element grid. */
template <typename ValueType_>
class Function
{
public:
    typedef ValueType_ ValueType;
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

    virtual ~Function() {}

    /** \brief Dimension of the space containing the grid on which the function
     *  is defined. */
    virtual int worldDimension() const = 0;

    /** \brief Number of components of the function's values.
     *
     *  E.g. 1 for a scalar-valued function. */
    virtual int codomainDimension() const = 0;

    /** \brief Retrieve types of geometrical data on which the function depend.
     *
     *  An implementation of this method for a particular Function subclass
     *  should modify the \p geomDeps bitfield by adding to it, using the
     *  bitwise OR operation, an appropriate combination of the flags defined
     *  in the enum GeometricalDataType.
     *
     *  For example, a subclass representing a function whose value at a point
     *  depends on its global coordinates and on the orientation of the vector
     *  normal to the element at that point should modify the argument as
     *  follows:
        \code
        geomDeps |= GLOBALS | NORMALS;
        \endcode
     */
    virtual void addGeometricalDependencies(size_t& geomDeps) const = 0;

    /** \brief Evaluate the function at a list of points.
     *
     *  \param[in] geomData
     *    Geometrical data related to \f$n \geq 0\f$ points on an element of a
     *    grid. The number of points, \f$n\f$, can be obtained by calling
     *    <tt>geomData.pointCount()</tt>.
     *  \param[out] result
     *    A 2-dimensional array intended to store the function
     *    values. On output, <tt>result(i, j)</tt> should contain the
     *    <em>i</em>th component of the function at the <em>j</em>th point.
     *
     *  An implementation of this method may assume that \p geomData contain
     *  all the types of geometrical data specified in the implementation of
     *  addGeometricalDependencies(). Before filling the array \p result, this
     *  method must ensure that its size is correct, calling
     *  <tt>arma::Mat::set_size()</tt> if necessary.
     */
    virtual void evaluate(const GeometricalData<CoordinateType>& geomData,
                          arma::Mat<ValueType>& result) const = 0;
};

} // namespace Fiber

#endif
