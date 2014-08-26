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

#ifndef fiber_collection_of_kernels_hpp
#define fiber_collection_of_kernels_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"

#include <utility>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename T> class CollectionOf3dArrays;
template <typename T> class CollectionOf4dArrays;
template <typename CoordinateType> class GeometricalData;
/** \endcond */

/** \ingroup weak_form_elements
 *  \brief Collection of kernels.
 *
 *  This class represents a collection of one or more integral operator
 *  *kernels*. A kernel is a function mapping a pair of points \f$(x, y)\f$
 *  located on two, possibly identical, elements of a grid to a scalar, vector
 *  or tensor of a fixed dimension, with real or complex elements. It can
 *  depend on any geometrical data related to \f$x\f$ and \f$y\f$---not only
 *  their global coordinates, but also the unit vectors normal to the grid at
 *  these points or the Jacobian matrices. All kernels in a collection are
 *  evaluated together and hence may reuse results of any intermediate
 *  calculations.
 *
 *  \tparam ValueType_
 *    Type used to represent the values of kernels (or the values of their
 *    components, in case or vector- or tensor-valued kernels). May be one of:
 *    <tt>float</tt>, <tt>double</tt>, <tt>std::complex<float></tt> or
 *    <tt>std::complex<double></tt>.
 */
template <typename ValueType_> class CollectionOfKernels {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  /** \brief Destructor */
  virtual ~CollectionOfKernels() {}

  /** \brief Retrieve types of geometrical data on which the kernels depend.
   *
   *  An implementation of this function for a particular kernel collection
   *  should modify the \p testGeomDeps and \p trialGeomDeps bitfields by
   *  adding to them, using the bitwise OR operation, an appropriate
   *  combination of the flags defined in the enum GeometricalDataType.
   *
   *  For example, a collection of kernels depending on the global
   *  coordinates of test and trial points and on the orientation of the
   *  vector normal to the trial element at trial points should modify the
   *  arguments as follows:
      \code
      testGeomDeps |= GLOBALS;
      trialGeomDeps |= GLOBALS | NORMALS;
      \endcode */
  virtual void addGeometricalDependencies(size_t &testGeomDeps,
                                          size_t &trialGeomDeps) const = 0;

  /** \brief Evaluate the kernels at a list of (test point, trial point) pairs.
   *
   *  \param[in] testGeomData
   *    Geometrical data related to \f$n \geq 0\f$ points on the test element.
   *    The number of points, \f$n\f$, can be obtained by calling
   *    <tt>testGeomData.pointCount()</tt>.
   *  \param[in] trialGeomData
   *    Geometrical data related to \f$n \geq 0\f$ points on the trial element.
   *    The number of test and trial points should be the same.
   *  \param[out] result
   *    A collection of 3-dimensional arrays intended to store the kernel
   *    values. On output, <tt>result[i][(j, k, p)</tt> should contain the
   *    (<em>j</em>, <em>k</em>)th entry in the tensor being the value of
   *    <em>i</em>th kernel at <em>p</em>th test point and <em>p</em>th trial
   *    point.
   *
   *  An implementation of this function may assume that \p testGeomData and
   *  \p trialGeomData contain all the types of geometrical data specified in
   *  the implementation of addGeometricalDependencies(). Before filling the
   *  arrays from \p result, this function must ensure that size of the array
   *  collection and the dimensions of the individual arrays are correct,
   *  calling <tt>CollectionOf3dArrays::set_size()</tt> and
   *  <tt>_3dArray::set_size()</tt> if necessary.
   */
  virtual void
  evaluateAtPointPairs(const GeometricalData<CoordinateType> &testGeomData,
                       const GeometricalData<CoordinateType> &trialGeomData,
                       CollectionOf3dArrays<ValueType> &result) const = 0;

  /** \brief Evaluate the kernels on a tensor grid of test and trial points.
   *
   *  \param[in] testGeomData
   *    Geometrical data related to \f$m \geq 0\f$ points on the test element.
   *    The number of points, \f$m\f$, can be obtained by calling
   *    <tt>testGeomData.pointCount()</tt>.
   *  \param[in] trialGeomData
   *    Geometrical data related to \f$n \geq 0\f$ points on the trial element.
   *    The number of points, \f$n\f$, can be obtained by calling
   *    <tt>trialGeomData.pointCount()</tt>.
   *  \param[out] result
   *    A collection of 4-dimensional arrays intended to store the kernel
   *    values. On output, <tt>result[i][(j, k, p, q)</tt> should contain the
   *    (<em>j</em>, <em>k</em>)th entry in the tensor being the value of
   *    <em>i</em>th kernel at <em>p</em>th test point and <em>q</em>th trial
   *    point.
   *
   *  An implementation of this function may assume that \p testGeomData and
   *  \p trialGeomData contain all the types of geometrical data specified in
   *  the implementation of addGeometricalDependencies(). Before filling the
   *  arrays from \p result, this function must ensure that size of the array
   *  collection and the dimensions of the individual arrays are correct,
   *  calling <tt>CollectionOf4dArrays::set_size()</tt> and
   *  <tt>_4dArray::set_size()</tt> if necessary.
   */
  virtual void
  evaluateOnGrid(const GeometricalData<CoordinateType> &testGeomData,
                 const GeometricalData<CoordinateType> &trialGeomData,
                 CollectionOf4dArrays<ValueType> &result) const = 0;

  /** \brief Currently unused. */
  virtual std::pair<const char *, int> evaluateClCode() const {
    throw std::runtime_error("CollectionOfKernels::evaluateClCode(): "
                             "not implemented");
  }

  virtual CoordinateType
  estimateRelativeScale(CoordinateType distance) const = 0;
};

} // namespace Fiber

#endif
