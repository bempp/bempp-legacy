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

#ifndef fiber_geometrical_data_hpp
#define fiber_geometrical_data_hpp

#include "../common/common.hpp"

#include "types.hpp"
#include "_3d_array.hpp"

#include <cassert>

namespace Fiber {

/** \brief Types of geometrical data. */
enum GeometricalDataType {
  GLOBALS = 0x0001,
  INTEGRATION_ELEMENTS = 0x0002,
  NORMALS = 0x0004,
  JACOBIANS_TRANSPOSED = 0x0008,
  JACOBIAN_INVERSES_TRANSPOSED = 0x0010,
  DOMAIN_INDEX = 0x0020
};

/** \cond FORWARD_DECL */
template <typename CoordinateType> class ConstGeometricalDataSlice;
/** \endcond */

/** \brief Storage of geometrical data.
 *
 *  \see Bempp::Geometry for a description of the data format (in particular,
 *  array ordering).
 */
template <typename CoordinateType> class GeometricalData {
public:
  Matrix<CoordinateType> globals;
  RowVector<CoordinateType> integrationElements;
  Fiber::_3dArray<CoordinateType> jacobiansTransposed;
  Fiber::_3dArray<CoordinateType> jacobianInversesTransposed;
  Matrix<CoordinateType> normals;
  int domainIndex;

  // For the time being, I (somewhat dangerously) assume that
  // integrationElements or globals or normals are always used
  int pointCount() const {
    int result = std::max(std::max(globals.cols(), normals.cols()),
                          integrationElements.cols());
    assert(result > 0);
    return result;
  }

  int dimWorld() const {
    const int dim = 3;
    return dim;
  }

  ConstGeometricalDataSlice<CoordinateType> const_slice(int point) const {
    return ConstGeometricalDataSlice<CoordinateType>(*this, point);
  }
};

/** \brief Access to slices of geometrical data.
 *
 *  This class gives access to a "slice" of geometrical data stored in a
 *  GeometricalData object, corresponding to a single point.
 */
template <typename CoordinateType> class ConstGeometricalDataSlice {
public:
  ConstGeometricalDataSlice(const GeometricalData<CoordinateType> &geomData,
                            int point)
      : m_geomData(geomData), m_point(point) {}

  CoordinateType global(int dim) const {
    return m_geomData.globals(dim, m_point);
  }
  CoordinateType integrationElement() const {
    return m_geomData.integrationElements(m_point);
  }
  CoordinateType jacobianTransposed(int row, int col) const {
    return m_geomData.jacobiansTransposed(row, col, m_point);
  }
  CoordinateType jacobianInverseTransposed(int row, int col) const {
    return m_geomData.jacobianInversesTransposed(row, col, m_point);
  }
  CoordinateType normal(int dim) const {
    return m_geomData.normals(dim, m_point);
  }
  int domainIndex() const { return m_geomData.domainIndex; }
  int dimWorld() const { return m_geomData.dimWorld(); }

  // Inefficient, but safe
  GeometricalData<CoordinateType> asGeometricalData() const {
    GeometricalData<CoordinateType> result;
    if (!is_empty(m_geomData.globals))
      result.globals = m_geomData.globals.col(m_point);
    if (!is_empty(m_geomData.integrationElements)) {
      result.integrationElements.resize(1);
      result.integrationElements(0) = m_geomData.integrationElements(m_point);
    }
    if (!m_geomData.jacobiansTransposed.is_empty()) {
      result.jacobiansTransposed.set_size(
          1, 1, m_geomData.jacobiansTransposed.extent(2));
      for (size_t i = 0; i < result.jacobiansTransposed.extent(2); ++i)
        result.jacobiansTransposed(0, 0, i) =
            m_geomData.jacobiansTransposed(m_point, m_point, i);
    }
    if (!m_geomData.jacobianInversesTransposed.is_empty()) {
      result.jacobianInversesTransposed.set_size(
          1, 1, m_geomData.jacobianInversesTransposed.extent(2));
      for (size_t i = 0; i < result.jacobianInversesTransposed.extent(2); ++i)
        result.jacobianInversesTransposed(0, 0, i) =
            m_geomData.jacobianInversesTransposed(m_point, m_point, i);
    }
    if (!is_empty(m_geomData.normals))
      result.normals = m_geomData.normals.col(m_point);
    result.domainIndex = m_geomData.domainIndex;
    return result;
  }

private:
  const GeometricalData<CoordinateType> &m_geomData;
  int m_point;
};

} // namespace Fiber

#endif
