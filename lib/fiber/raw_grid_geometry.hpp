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

#ifndef fiber_raw_grid_geometry_hpp
#define fiber_raw_grid_geometry_hpp

#include "../common/common.hpp"

#include "types.hpp"

namespace Fiber {

template <typename CoordinateType> class RawGridGeometry {
public:
  RawGridGeometry(int gridDim, int worldDim)
      : m_gridDim(gridDim), m_worldDim(worldDim) {
    if (gridDim > worldDim)
      throw std::invalid_argument("RawGridGeometry::RawGridGeometry(): "
                                  "grid dimension cannot be larger than "
                                  "world dimension");
  }

  // Const accessors

  const Matrix<CoordinateType> &vertices() const { return m_vertices; }

  const Matrix<int> &elementCornerIndices() const {
    return m_elementCornerIndices;
  }

  const Matrix<char> &auxData() const { return m_auxData; }

  const std::vector<int> &domainIndices() const { return m_domainIndices; }

  int elementCount() const { return m_elementCornerIndices.cols(); }

  int gridDimension() const { return m_gridDim; }

  int worldDimension() const { return m_worldDim; }

  /** \brief Indices of the corners of the given element. */
  Vector<int> elementCornerIndices(int elementIndex) const {
    const int n = elementCornerCount(elementIndex);
    // return m_elementCornerIndices(arma::span(0, n - 1),
    //                              arma::span(elementIndex));
    return m_elementCornerIndices.block(0, elementIndex, n, 1);
  }

  /** \brief Number of corners of the given element. */
  int elementCornerCount(int elementIndex) const {
    int n = m_elementCornerIndices.rows();
    while (m_elementCornerIndices(n - 1, elementIndex) < 0)
      --n;
    return n;
  }

  /** \brief Domain index of the given element. */
  int domainIndex(int elementIndex) const {
    return m_domainIndices[elementIndex];
  }

  // Non-const accessors (currently needed for construction)

  Matrix<CoordinateType> &vertices() { return m_vertices; }

  Matrix<int> &elementCornerIndices() { return m_elementCornerIndices; }

  Matrix<char> &auxData() { return m_auxData; }

  std::vector<int> &domainIndices() { return m_domainIndices; }

  // Auxiliary functions

  template <typename Geometry>
  void setupGeometry(int elementIndex, Geometry &geometry) const {
    const int dimGrid = m_vertices.rows();
    size_t cornerCount = 0;
    for (; cornerCount < m_elementCornerIndices.rows(); ++cornerCount)
      if (m_elementCornerIndices(cornerCount, elementIndex) < 0)
        break;
    Matrix<CoordinateType> corners(dimGrid, cornerCount);
    for (size_t cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex)
      corners.col(cornerIndex) =
          m_vertices.col(m_elementCornerIndices(cornerIndex, elementIndex));
    Eigen::Map<Matrix<char>> auxDataMap(
        const_cast<char *>(m_auxData.col(elementIndex).data()),
        m_auxData.rows(), 1);
    geometry.setup(corners, auxDataMap);
  }

private:
  int m_gridDim;
  int m_worldDim;
  Matrix<CoordinateType> m_vertices;
  Matrix<int> m_elementCornerIndices;
  Matrix<char> m_auxData;
  std::vector<int> m_domainIndices;
};

} // namespace Fiber

#endif
