// Copyright (C) 2011 by the BEM++ Authors
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

#ifndef bempp_element_hpp
#define bempp_element_hpp

#include "utils.hpp"

/**
 * Abstract basis class for elements
 */

class Element
{
public:
  /** Constructor */
  Element() {}
  /** Destructor */
  virtual ~Element() {}
  /**
   * Calculate the Jacobian of the coordinate transform mapping
   * logical coordinates (xi, eta) into physical coordinates (x, y)
   * at specified points.
   *
   * @param[in]  n       number of points
   * @param[in]  pts     list of points (logical coordinates)
   * @param[out] result  Jacobian at the specified points
   */
  virtual void jacobian(int n, const Point* pts, double* result) const = 0;
};

/**
 * Triangular element
 */
class TriangularElement : public Element
{
public:
  virtual void jacobian(int n, const Point* pts, double* result) const;
  /** Triangle's vertices (physical coordinates) */
  Point vertices[3];
};

/**
 * Quadrilateral element
 */
class QuadrilateralElement : public Element
{
public:
  virtual void jacobian(int n, const Point* pts, double* result) const;
  /** Quadrilateral's vertices (physical coordinates) */
  Point vertices[4];
};

#endif
