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
#include "../space/space.hpp"

#include "../common/not_implemented_error.hpp"

namespace Bempp {

namespace {

double min3(double x, double y, double z) {
  return std::min(x, std::min(y, z));
}

double max3(double x, double y, double z) {
  return std::max(x, std::max(y, z));
}

bool isNew(Vector<double> &intersection,
           const std::vector<Vector<double>> &intersections) {
  const double EPSILON = 1e-10;
  for (size_t i = 0; i < intersections.size(); ++i)
    if (fabs(intersections[i](0) - intersection(0)) < EPSILON &&
        fabs(intersections[i](1) - intersection(1)) < EPSILON &&
        fabs(intersections[i](2) - intersection(2)) < EPSILON)
      return false;
  return true;
}

} // namespace

void Grid::signalGridUpdate() const
{

    gridUpdateSignal();

}

boost::signals2::connection Grid::connect(const std::function<void()>& f) const
{
    return gridUpdateSignal.connect(f);
}

bool Grid::isBarycentricRepresentationOf(const Grid &other) const {
  if (!other.hasBarycentricGrid())
    return false;
  else
    return (this == other.barycentricGrid().get());
}

void Grid::getBoundingBox(Vector<double> &lowerBound,
                          Vector<double> &upperBound) const {
  // In this simple implementation we assume that all elements are flat.
  if (m_lowerBound.rows() == dimWorld() && m_upperBound.rows() == dimWorld()) {
    lowerBound = m_lowerBound;
    upperBound = m_upperBound;
    return;
  }

  std::unique_ptr<GridView> view = leafView();

  Matrix<double> vertices;
  Matrix<int> elementCorners; // unused
  Matrix<char> auxData;       // unused
  view->getRawElementData(vertices, elementCorners, auxData);

  m_lowerBound = lowerBound =
      vertices.rowwise().minCoeff(); // min. value in each row
  m_upperBound = upperBound =
      vertices.rowwise().maxCoeff(); // max. value in each row
}

std::vector<bool> areInside(const Grid &grid, const Matrix<double> &points) {
  if (grid.dim() != 2 || grid.dimWorld() != 3)
    throw NotImplementedError("areInside(): currently implemented only for"
                              "2D grids embedded in 3D spaces");

  std::unique_ptr<GridView> view = grid.leafView();
  std::unique_ptr<EntityIterator<0>> it = view->entityIterator<0>();

  std::vector<Matrix<double>> triangles;
  triangles.reserve(view->entityCount(0));

  {
    Matrix<double> triangle(grid.dimWorld(), 3);
    Matrix<double> corners;
    while (!it->finished()) {
      const Entity<0> &element = it->entity();
      element.geometry().getCorners(corners);
      if (corners.cols() == 3) // triangle
        triangles.push_back(corners);
      else if (corners.cols() == 4) { // quadrilateral, split into 2 triangles
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
  RowVector<double> xMins(triangleCount);
  RowVector<double> xMaxs(triangleCount);
  RowVector<double> yMins(triangleCount);
  RowVector<double> yMaxs(triangleCount);
  RowVector<double> zMins(triangleCount);
  RowVector<double> zMaxs(triangleCount);
  for (size_t i = 0; i < triangleCount; ++i) {
    const Matrix<double> &triangle = triangles[i];
    xMins(i) = min3(triangle(0, 0), triangle(0, 1), triangle(0, 2));
    xMaxs(i) = max3(triangle(0, 0), triangle(0, 1), triangle(0, 2));
    yMins(i) = min3(triangle(1, 0), triangle(1, 1), triangle(1, 2));
    yMaxs(i) = max3(triangle(1, 0), triangle(1, 1), triangle(1, 2));
    zMins(i) = min3(triangle(2, 0), triangle(2, 1), triangle(2, 2));
    zMaxs(i) = max3(triangle(2, 0), triangle(2, 1), triangle(2, 2));
  }

  const double gridXMin = xMins.minCoeff();
  const double gridXMax = xMaxs.maxCoeff();
  const double gridYMin = yMins.minCoeff();
  const double gridYMax = yMaxs.maxCoeff();
  const double gridZMin = zMins.minCoeff();
  const double gridZMax = zMaxs.maxCoeff();

  const size_t pointCount = points.cols();
  std::vector<bool> result(pointCount, false);
  std::vector<Vector<double>> intersections;
  Vector<double> intersection(3);
  for (size_t pt = 0; pt < pointCount; ++pt) {
    if (points(0, pt) < gridXMin || points(0, pt) > gridXMax ||
        points(1, pt) < gridYMin || points(1, pt) > gridYMax ||
        points(2, pt) < gridZMin || points(2, pt) > gridZMax)
      continue; // point outside grid's bounding box
    intersections.clear();
    for (size_t tri = 0; tri < triangleCount; ++tri)
      if (points(0, pt) >= xMins(tri) && points(0, pt) <= xMaxs(tri) &&
          points(1, pt) >= yMins(tri) && points(1, pt) <= yMaxs(tri))
        if (zRayIntersectsTriangle(
                points.col(pt).data(), triangles[tri].col(0).data(),
                triangles[tri].col(1).data(), triangles[tri].col(2).data(),
                intersection.col(0).data()) > 0.)
          if (isNew(intersection, intersections))
            intersections.push_back(intersection);
    result[pt] = intersections.size();
  }
  return result;
}

std::vector<bool> areInside(const Grid &grid, const Matrix<float> &points) {
  Matrix<double> pointsDouble(points.rows(), points.cols());
  for (int c = 0; c < points.cols(); ++c)
    for (int r = 0; r < points.rows(); ++r)
      pointsDouble(r, c) = points(r, c);
  return areInside(grid, pointsDouble);
}

unsigned int Grid::elementInsertionIndex(const Entity<0> &element) const {
  throw std::runtime_error("Grid::elementInsertionIndex(): "
                           "method not implemented.");
}

unsigned int Grid::vertexInsertionIndex(const Entity<2> &vertex) const {
  throw std::runtime_error("Grid::vertexInsertionIndex(): "
                           "method not implemented.");
}

} // namespace Bempp
