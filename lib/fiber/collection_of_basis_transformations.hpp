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

#ifndef fiber_collection_of_basis_transformations_hpp
#define fiber_collection_of_basis_transformations_hpp

#include "../common/common.hpp"

#include "scalar_traits.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename T> class CollectionOf3dArrays;
template <typename ValueType> class BasisData;
template <typename CoordinateType> class GeometricalData;
/** \endcond */

/** \brief Collection of kernels.
 *
 *  This class represents a collection of one or more basis function
 *  transformations. A basis function transformation is a function mapping a
 *  point \f$x\f$ located at an element of a grid to a scalar or vector
 *  of a fixed dimension, with real or complex elements. In addition to any
 *  geometrical data related to \f$x\f$, such as its global coordinates or the
 *  unit vectors normal to the grid at \f$x\f$, it can depend on the value
 *  and/or the first derivatives of a basis function defined on the reference
 *  element (e.g.\ the unit triangle or the unit square). Example basis
 *  function transformations include the mapping of basis functions to shape
 *  functions---expressed with the identity mapping for scalar basis functions
 *  or the Piola transform for vector basis function---and the mapping of basis
 *  functions to the surface curls of shape functions. All basis function
 *  transformations in a collection are evaluated together and hence may reuse
 *  results of any intermediate calculations.
 *
 *  A basis function transformation is assumed to preserve the type of the
 *  values of basis functions, i.e. it should map real-valued basis functions
 *  to real-valued functions and complex-valued basis functions to
 *  complex-valued functions.
 *
 *  \tparam CoordinateType_
 *    Type used to represent coordinates (either <tt>float</tt> or
 *    <tt>double</tt>). */
template <typename CoordinateType_>
class CollectionOfBasisTransformations
{
public:
    typedef CoordinateType_ CoordinateType;
    typedef typename ScalarTraits<CoordinateType>::ComplexType ComplexType;

    /** \brief Destructor. */
    virtual ~CollectionOfBasisTransformations()
    {}

    /** \brief Return the number of transformations belonging to the collection. */
    virtual int transformationCount() const = 0;

    /** \brief Return the number of components of basis functions acted upon.
     *
     *  For instance, the implementation of this function for a collection of
     *  transformations operating on scalar basis functions should return 1; for
     *  a collection operating on vector-valued basis functions with three
     *  components, this function should return 3. */
    virtual int argumentDimension() const = 0;

    /** \brief Return the number of components of the result of i'th
     *  transformation.
     *
     *  Transformation indices start from 0.
     *
     *  For example, if the first transformation produces a scalar function,
     *  the implementation of this function should return 1 if \p i == 0.
     *
     *  The behaviour of this function for \p i < 0 or \p >= transformationCount()
     *  is not specified. */
    virtual int resultDimension(int i) const = 0;

    /** \brief Retrieve the types of data on which the transformations depend.
     *
     *  An implementation of this function should modify the \p basisDeps and
     *  \p geomDeps bitfields by adding to them, using the bitwise OR
     *  operation, an appropriate combination of the flags defined in the enums
     *  BasisDataType and GeometricalDataType.
     *
     *  For example, a collection of transformations depending on the values
     *  and first derivatives of basis functions and the global coordinates of
     *  points at which the transformed functions are evaluated should modify
     *  the arguments as follows:
     *
        \code
        basisDeps |= VALUES | DERIVATIVES;
        geomDeps |= GLOBALS;
        \endcode */
    virtual void addDependencies(size_t& basisDeps, size_t& geomDeps) const = 0;

    /** \brief Evaluate transformations of real-valued basis functions.
     *
     *  \param[in] basisData
     *    Values and/or derivatives of \f$m \geq 0 \f$ basis functions at \f$n
     *    \geq 0\f$ points on an element. The number of points, \f$n\f$, can be
     *    obtained by calling <tt>basisData.pointCount()</tt>; the number of
     *    basis functions, \f$m\f$, can be obtained by calling
     *    <tt>basisData.functionCount()</tt>.
     *  \param[in] geomData
     *    Geometrical data related to the \f$n\f$ points at which the
     *    basis functions have been evaluated.
     *  \param[out] result
     *    A collection of 3-dimensional arrays intended to store the
     *    transformed basis function values. On output, <tt>result[i][(j, k,
     *    p)</tt> should contain the <em>j</em> element of the vector being the
     *    value of <em>i</em>th transformation of <em>k</em>th basis function
     *    at <em>p</em>th point.
     *
     *  An implementation of this function may assume that \p basisData and \p
     *  geomData contain all the types of data specified in the implementation
     *  of addDependencies(). Before filling the arrays from \p result, this
     *  function must ensure that size of the array collection and the
     *  dimensions of the individual arrays are correct, calling
     *  <tt>CollectionOf3dArrays::set_size()</tt> and
     *  <tt>_3dArray::set_size()</tt> if necessary.
     */
    void evaluate(
            const BasisData<CoordinateType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<CoordinateType>& result) const {
        evaluateImplReal(basisData, geomData, result);
    }

    /** \brief Evaluate transformations of complex-valued basis functions.
     *
     *  See the documentation of the other overload for the description of
     *  function parameters.
     */
    void evaluate(
            const BasisData<ComplexType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<ComplexType>& result) const {
        evaluateImplComplex(basisData, geomData, result);
    }

private:
    /** \brief Evaluate transformations of real-valued basis functions.
     *
     *  \param[in] basisData
     *    Values and/or derivatives of \f$m \geq 0 \f$ basis functions at \f$n
     *    \geq 0\f$ points on an element. The number of points, \f$n\f$, can be
     *    obtained by calling <tt>basisData.pointCount()</tt>; the number of
     *    basis functions, \f$m\f$, can be obtained by calling
     *    <tt>basisData.functionCount()</tt>.
     *  \param[in] geomData
     *    Geometrical data related to the \f$n\f$ points at which the
     *    basis functions have been evaluated.
     *  \param[out] result
     *    A collection of 3-dimensional arrays intended to store the
     *    transformed basis function values. On output, <tt>result[i][(j, k,
     *    p)</tt> should contain the <em>j</em> element of the vector being the
     *    value of <em>i</em>th transformation of <em>k</em>th basis function
     *    at <em>p</em>th point.
     *
     *  This is a pure virtual function that must be overridden in subclasses
     *  of CollectionOfBasisTransformations.
     *
     *  An implementation of this function may assume that \p basisData and \p
     *  geomData contain all the types of data specified in the implementation
     *  of addDependencies(). Before filling the arrays from \p result, this
     *  function must ensure that size of the array collection and the
     *  dimensions of the individual arrays are correct, calling
     *  <tt>CollectionOf3dArrays::set_size()</tt> and
     *  <tt>_3dArray::set_size()</tt> if necessary.
     */
    virtual void evaluateImplReal(
            const BasisData<CoordinateType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<CoordinateType>& result) const = 0;

    /** \brief Evaluate transformations of complex-valued basis functions.
     *
     *  See the documentation of evaluateImplReal() for the description of
     *  function parameters. */
    virtual void evaluateImplComplex(
            const BasisData<ComplexType>& basisData,
            const GeometricalData<CoordinateType>& geomData,
            CollectionOf3dArrays<ComplexType>& result) const = 0;
};

} // namespace Fiber

#endif

