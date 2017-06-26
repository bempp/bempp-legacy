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

#ifndef bempp_bounding_box_helpers_hpp
#define bempp_bounding_box_helpers_hpp

#include "../common/eigen_support.hpp"
#include "bounding_box.hpp"
#include "common.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace Bempp {

/** \brief Create bounding box with reference set to mid point */
template <typename CoordinateType>
BoundingBox<CoordinateType> createBoundingBox(
    CoordinateType xmin, CoordinateType ymin, CoordinateType zmin,
    CoordinateType xmax, CoordinateType ymax, CoordinateType zmax)
{

    BoundingBox<CoordinateType> bbox;
    bbox.lbound.x = xmin;
    bbox.lbound.y = ymin;
    bbox.lbound.z = zmin;
    bbox.ubound.x = xmax;
    bbox.ubound.y = ymax;
    bbox.ubound.z = zmax;
    bbox.reference.x = (xmax - xmin) / 2;
    bbox.reference.y = (ymax - ymin) / 2;
    bbox.reference.z = (zmax - zmin) / 2;
    return bbox;
}

/** \relates BoundingBox
 *  \brief Extend the bounding box \p bbox to include all points contained in
 *  the columns of \p points.
 *
 *  \note The array \p points needs to have three rows.
 */
template <typename CoordinateType>
void extendBoundingBox(BoundingBox<CoordinateType>& bbox,
    const Matrix<CoordinateType>& points)
{
    assert(points.rows() == 3);
    for (size_t j = 0; j < points.cols(); ++j) {
        bbox.lbound.x = std::min(bbox.lbound.x, points(0, j));
        bbox.lbound.y = std::min(bbox.lbound.y, points(1, j));
        bbox.lbound.z = std::min(bbox.lbound.z, points(2, j));
        bbox.ubound.x = std::max(bbox.ubound.x, points(0, j));
        bbox.ubound.y = std::max(bbox.ubound.y, points(1, j));
        bbox.ubound.z = std::max(bbox.ubound.z, points(2, j));
    }
}

/** \relates BoundingBox
 *  \brief Set the reference point of the bounding box \p bbox to \p point.
 *
 *  \note \p point must be a three-component column vector.
 */
template <typename CoordinateType>
void setBoundingBoxReference(BoundingBox<CoordinateType>& bbox,
    const Vector<CoordinateType>& point)
{
    assert(point.rows() == 3);
    bbox.reference.x = point(0);
    bbox.reference.y = point(1);
    bbox.reference.z = point(2);
}

template <typename CoordinateType>
Vector<CoordinateType> getBoundingBoxSize(const BoundingBox<CoordinateType>& boundingBox)
{
    Vector<CoordinateType> result;
    result[0] = boundingBox.ubound.x - boundingBox.lbound.x;
    result[1] = boundingBox.ubound.y - boundingBox.lbound.y;
    result[2] = boundingBox.ubound.z - boundingBox.lbound.z;

    return result;
}

} // namespace Bempp

#endif
