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

#include "common.hpp"
#include "armadillo_fwd.hpp"
#include "bounding_box.hpp"

#include <cassert>
#include <cmath>

namespace Bempp
{

/** \relates BoundingBox
 *  \brief Extend the bounding box \p bbox to include all points contained in
 *  the columns of \p points.
 *
 *  \note The array \p points needs to have three rows.
 */
template <typename CoordinateType>
void extendBoundingBox(BoundingBox<CoordinateType>& bbox,
                       const arma::Mat<CoordinateType>& points)
{
    assert(points.n_rows == 3);
    for (size_t j = 0; j < points.n_cols; ++j) {
        bbox.lbound.x =
            std::min(bbox.lbound.x, points(0, j));
        bbox.lbound.y =
            std::min(bbox.lbound.y, points(1, j));
        bbox.lbound.z =
            std::min(bbox.lbound.z, points(2, j));
        bbox.ubound.x =
            std::max(bbox.ubound.x, points(0, j));
        bbox.ubound.y =
            std::max(bbox.ubound.y, points(1, j));
        bbox.ubound.z =
            std::max(bbox.ubound.z, points(2, j));
    }
}

/** \relates BoundingBox
 *  \brief Set the reference point of the bounding box \p bbox to \p point.
 *
 *  \note \p point must be a three-component column vector.
 */
template <typename CoordinateType>
void setBoundingBoxReference(BoundingBox<CoordinateType>& bbox,
                             const arma::Col<CoordinateType>& point)
{
    assert(point.n_rows == 3);
    bbox.reference.x = point(0);
    bbox.reference.y = point(1);
    bbox.reference.z = point(2);
}

} // namespace Bempp

#endif
