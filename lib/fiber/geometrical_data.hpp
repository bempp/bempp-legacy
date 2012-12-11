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

#include "../common/armadillo_fwd.hpp"

#include <cassert>

namespace Fiber
{

/** \brief Types of geometrical data. */
enum GeometricalDataType
{
    GLOBALS = 0x0001,
    INTEGRATION_ELEMENTS = 0x0002,
    NORMALS = 0x0004,
    JACOBIANS_TRANSPOSED = 0x0008,
    JACOBIAN_INVERSES_TRANSPOSED = 0x0010
};

/** \cond FORWARD_DECL */
template <typename CoordinateType> class ConstGeometricalDataSlice;
/** \endcond */

template <typename CoordinateType>
struct GeometricalData
{
    arma::Mat<CoordinateType> globals;
    arma::Row<CoordinateType> integrationElements;
    arma::Cube<CoordinateType> jacobiansTransposed;
    arma::Cube<CoordinateType> jacobianInversesTransposed;
    arma::Mat<CoordinateType> normals;

    // For the time being, I (somewhat dangerously) assume that
    // integrationElements or globals or normals are always used
    int pointCount() const {
        int result = std::max(std::max(globals.n_cols, normals.n_cols),
                                       integrationElements.n_cols);
        assert(result > 0);
        return result;
    }

    int dimWorld() const {
        int result = std::max(globals.n_rows, normals.n_rows);
        assert(result > 0);
        return result;
    }

    ConstGeometricalDataSlice<CoordinateType> const_slice(int point) const {
        return ConstGeometricalDataSlice<CoordinateType>(*this, point);
    }
};

template <typename CoordinateType>
class ConstGeometricalDataSlice
{
public:
    ConstGeometricalDataSlice(const GeometricalData<CoordinateType>& geomData,
                              int point) :
        m_geomData(geomData), m_point(point) {}

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

    int dimWorld() const {
        return m_geomData.dimWorld();
    }

    // Inefficient, but safe
    GeometricalData<CoordinateType> asGeometricalData() const {
        GeometricalData<CoordinateType> result;
        if (!m_geomData.globals.is_empty())
            result.globals = m_geomData.globals.col(m_point);
        if (!m_geomData.integrationElements.is_empty())
            result.integrationElements =
                    m_geomData.integrationElements(m_point);
        if (!m_geomData.jacobiansTransposed.is_empty())
            result.jacobiansTransposed =
                    m_geomData.jacobiansTransposed.slices(m_point, m_point);
        if (!m_geomData.jacobianInversesTransposed.is_empty())
            result.jacobianInversesTransposed =
                    m_geomData.jacobianInversesTransposed.slices(m_point, m_point);
        if (!m_geomData.normals.is_empty())
            result.normals = m_geomData.normals.col(m_point);
        return result;
    }

private:
    const GeometricalData<CoordinateType>& m_geomData;
    int m_point;
};

} // namespace Fiber

#endif
