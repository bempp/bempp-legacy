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

#include "grid.hpp"

#include "entity.hpp"
#include "entity_iterator.hpp"
#include "geometry.hpp"
#include "grid_view.hpp"
#include "ray_triangle_intersection.hpp"

#include "../common/not_implemented_error.hpp"

namespace Bempp
{

namespace {

double min3(double x, double y, double z)
{
    return std::min(x, std::min(y, z));
}

double max3(double x, double y, double z)
{
    return std::max(x, std::max(y, z));
}

bool isNew(const arma::Col<double>& intersection,
           const std::vector<arma::Col<double> >& intersections)
{
    const double EPSILON = 1e-10;
    for (size_t i = 0; i < intersections.size(); ++i)
        if (fabs(intersections[i](0) - intersection(0)) < EPSILON &&
                fabs(intersections[i](1) - intersection(1)) < EPSILON &&
                fabs(intersections[i](2) - intersection(2)) < EPSILON)
            return false;
    return true;
}

} // namespace

std::vector<bool> areInside(const Grid& grid, const arma::Mat<double>& points)
{
    if (grid.dim() != 2 || grid.dimWorld() != 3)
        throw NotImplementedError("areInside(): currently implemented only for"
                                  "2D grids embedded in 3D spaces");

    std::auto_ptr<GridView> view = grid.leafView();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();

    std::vector<arma::Mat<double> > triangles;
    triangles.reserve(view->entityCount(0));

    {
        arma::Mat<double> triangle(grid.dimWorld(), 3);
        arma::Mat<double> corners;
        while (!it->finished()) {
            const Entity<0>& element = it->entity();
            element.geometry().getCorners(corners);
            if (corners.n_cols == 3) // triangle
                triangles.push_back(corners);
            else if (corners.n_cols == 4) { // quadrilateral, split into 2 triangles
                // NOTE: this won't work for concave quads,
                triangle.col(0) = corners.col(0);
                triangle.col(1) = corners.col(1);
                triangle.col(2) = corners.col(2);
                triangles.push_back(triangle);
                triangle.col(0) = corners.col(2);
                triangle.col(1) = corners.col(3);
                triangle.col(2) = corners.col(0);
                triangles.push_back(triangle);
            } else
                throw std::runtime_error("areInside(): unknown element type");
            it->next();
        }
    }

    const size_t triangleCount = triangles.size();
    arma::Row<double> xMins(triangleCount);
    arma::Row<double> xMaxs(triangleCount);
    arma::Row<double> yMins(triangleCount);
    arma::Row<double> yMaxs(triangleCount);
    arma::Row<double> zMins(triangleCount);
    arma::Row<double> zMaxs(triangleCount);
    for (size_t i = 0; i < triangleCount; ++i) {
        const arma::Mat<double>& triangle = triangles[i];
        xMins(i) = min3(triangle(0, 0), triangle(0, 1), triangle(0, 2));
        xMaxs(i) = max3(triangle(0, 0), triangle(0, 1), triangle(0, 2));
        yMins(i) = min3(triangle(1, 0), triangle(1, 1), triangle(1, 2));
        yMaxs(i) = max3(triangle(1, 0), triangle(1, 1), triangle(1, 2));
        zMins(i) = min3(triangle(2, 0), triangle(2, 1), triangle(2, 2));
        zMaxs(i) = max3(triangle(2, 0), triangle(2, 1), triangle(2, 2));
    }

    const double gridXMin = arma::min(xMins);
    const double gridXMax = arma::max(xMaxs);
    const double gridYMin = arma::min(yMins);
    const double gridYMax = arma::max(yMaxs);
    const double gridZMin = arma::min(zMins);
    const double gridZMax = arma::max(zMaxs);

    const size_t pointCount = points.n_cols;
    std::vector<bool> result(pointCount, false);
    std::vector<arma::Col<double> > intersections;
    arma::Col<double> intersection(3);
    for (size_t pt = 0; pt < pointCount; ++pt) {
        if (points(0, pt) < gridXMin || points(0, pt) > gridXMax ||
                points(1, pt) < gridYMin || points(1, pt) > gridYMax ||
                points(2, pt) < gridZMin || points(2, pt) > gridZMax)
            continue; // point outside grid's bounding box
        intersections.clear();
        for (size_t tri = 0; tri < triangleCount; ++tri)
            if (points(0, pt) >= xMins(tri) && points(0, pt) <= xMaxs(tri) &&
                    points(1, pt) >= yMins(tri) && points(1, pt) <= yMaxs(tri))
                if (zRayIntersectsTriangle(points.colptr(pt),
                                           triangles[tri].colptr(0),
                                           triangles[tri].colptr(1),
                                           triangles[tri].colptr(2),
                                           intersection.colptr(0)) > 0.)
                    if (isNew(intersection, intersections))
                        intersections.push_back(intersection);
        result[pt] = intersections.size();
    }
    return result;
}

std::vector<bool> areInside(const Grid& grid, const arma::Mat<float>& points)
{
    arma::Mat<double> pointsDouble(points.n_rows, points.n_cols);
    for (int c = 0; c < points.n_cols; ++c)
        for (int r = 0; r < points.n_rows; ++r)
            pointsDouble(r, c) = points(r, c);
    return areInside(grid, pointsDouble);
}

} // namespace Bempp
