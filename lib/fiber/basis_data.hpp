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

#ifndef fiber_basis_data_hpp
#define fiber_basis_data_hpp

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"

#include "../common/multidimensional_arrays.hpp"

namespace Fiber {

enum BasisDataType { VALUES = 0x0001, DERIVATIVES = 0x0002 };

/** \cond FORWARD_DECL */
template <typename ValueType> class ConstBasisDataSlice;
/** \endcond */

/** \brief Storage of values and/or derivatives of shape functions. */
template <typename ValueType> class BasisData {
public:
  /** \brief Values of shape functions.
   *
   *  Data format:
   *  <tt>values(i, k, l) =</tt> \f$(f_k)_i(x_l)\f$, i.e.
   *  ith component of kth shape function at lth point. */
  _3dArray<ValueType> values;

  /** \brief Derivatives of shape functions.
   *
   *  Data format:
   *  <tt>derivatives(i, j, k, l) =</tt> \f$(d_j (f_k)_i)(x_l)\f$, i.e.
   *  the derivative in direction j of ith component of kth shape function at
   *  lth point. */
  _4dArray<ValueType> derivatives;

  /** \brief Return number of components of shape functions. */
  int componentCount() const {
    return std::max<int>(values.extent(0), derivatives.extent(0));
  }

  /** \brief Return number of shape functions. */
  int functionCount() const {
    return std::max<int>(values.extent(1), derivatives.extent(2));
  }

  /** \brief Return number of points at which the shape functions have been
   *  calculated. */
  int pointCount() const {
    return std::max<int>(values.extent(2), derivatives.extent(3));
  }

  /** \brief Return a constant slice of the data, corresponding to a given
   *  shape function and point. */
  ConstBasisDataSlice<ValueType> const_slice(int function, int point) const {
    return ConstBasisDataSlice<ValueType>(*this, function, point);
  }
};

/** \brief Access to values and/or derivatives of shape functions.
 *
 *  This class gives access to a "slice" of data stored in a BasisData object,
 *  corresponding to a single shape function and point.
 */
template <typename ValueType> class ConstBasisDataSlice {
public:
  /** \brief Constructor. */
  ConstBasisDataSlice(const BasisData<ValueType> &basisData, int function,
                      int point)
      : m_basisData(basisData), m_function(function), m_point(point) {}

  /** \brief Return the value of the \p dim'th component of the function. */
  ValueType values(int dim) const {
    return m_basisData.values(dim, m_function, m_point);
  }

  /** \brief Return the value of the \p dim'th component of the derivative of
   *  the function in a given direction. */
  ValueType derivatives(int dim, int direction) const {
    return m_basisData.derivatives(dim, direction, m_function, m_point);
  }

  /** \brief Return the number of components of the function. */
  int componentCount() const { return m_basisData.componentCount(); }

private:
  const BasisData<ValueType> &m_basisData;
  int m_function, m_point;
};

} // namespace Fiber

#endif // BASIS_DATA_TYPES_HPP
