// Copyright (C) 2011-2012 by the Fiber Authors
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

#include <armadillo>

namespace Fiber
{

template <typename ValueType>
class RawGridGeometry
{
public:
    // Const accessors

    const arma::Mat<ValueType>& vertices() const {
        return m_vertices;
    }

    const arma::Mat<int>& elementCornerIndices() const {
        return m_elementCornerIndices;
    }

    const arma::Mat<char>& auxData() const {
        return m_auxData;
    }

    /** \brief Indices of the corners of the given element. */
    arma::Col<int> elementCornerIndices(int elementIndex) const {
        const int n = elementCornerCount(elementIndex);
        return m_elementCornerIndices(arma::span(0, n - 1), arma::span(elementIndex));
    }

    /** \brief Number of corners of the given element. */
    int elementCornerCount(int elementIndex) const {
        int n = m_elementCornerIndices.n_rows;
        while (m_elementCornerIndices(n - 1, elementIndex) < 0)
            --n;
        return n;
    }

    // Non-const accessors (currently needed for construction)

    arma::Mat<ValueType>& vertices() {
        return m_vertices;
    }

    arma::Mat<int>& elementCornerIndices() {
        return m_elementCornerIndices;
    }

    arma::Mat<char>& auxData() {
        return m_auxData;
    }

    // Auxiliary functions

    template <typename Geometry>
    void setupGeometry(int elementIndex, Geometry& geometry) const
    {
        const int dimGrid = m_vertices.n_rows;
        int cornerCount = 0;
        for (; cornerCount < m_elementCornerIndices.n_rows; ++cornerCount)
            if (m_elementCornerIndices(cornerCount, elementIndex) < 0)
                break;
        arma::Mat<ValueType> corners(dimGrid, cornerCount);
        for (int cornerIndex = 0; cornerIndex < cornerCount; ++cornerIndex)
            corners.col(cornerIndex) = m_vertices.col(
                        m_elementCornerIndices(cornerIndex, elementIndex));
        geometry.setup(corners, m_auxData.unsafe_col(elementIndex));
    }

private:
    arma::Mat<ValueType> m_vertices;
    arma::Mat<int> m_elementCornerIndices;
    arma::Mat<char> m_auxData;
};

} // namespace Fiber

#endif
