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

// Code adapted from
// http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/

namespace Bempp
{

namespace {

/* dest = a - b */
inline void substract(double* dest, const double* a, const double* b)
{
    dest[0] = a[0] - b[0];
    dest[1] = a[1] - b[1];
    dest[2] = a[2] - b[2];
}

inline void crossProduct(double* dest, const double* v1, const double* v2)
{
    dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
    dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
    dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

inline double innerProduct(const double* v1, const double* v2)
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

} // namespace

// Tests if the ray (p[0], p[1], p[2]) + alpha (0, 0, 1) (alpha >= 0)
// intersects the triangle with vertices v0, v1 and v2. Returns 0 if
// it doesn't intersect, 1./6. if it intersects a vertex, 0.5 if it
// intersects an edge, and 1 if it passes through the interior of the
// triangle. If there is an intersection, its position is stored in the
// argument 'intersection'.
double zRayIntersectsTriangle(
        const double *p, const double *v0, const double *v1, const double *v2,
        double* intersection)
{
    double d[3] = {0., 0., 1.}; // ray directions
    double e1[3],e2[3],h[3],s[3],q[3];
    double a,f,t,u,v;
    substract(e1,v1,v0);
    substract(e2,v2,v0);

    crossProduct(h,d,e2);
    a = innerProduct(e1,h);

    const double eps = 1e-10;
    if (a > -eps && a < eps)
        return 0.;

    f = 1/a;
    substract(s,p,v0);
    u = f * (innerProduct(s,h));

    if (u < 0.0 || u > 1.0)
        return 0.;

    crossProduct(q,s,e1);
    v = f * innerProduct(d,q);

    if (v < 0.0 || u + v > 1.0)
        return 0.;

    // at this stage we can compute t to find out where
    // the intersection point is on the line
    t = f * innerProduct(e2,q);

    if (t < eps) // this means that there is a line intersection
                 // but not a ray intersection
        return 0.;
    else { // ray intersection
        for (int i = 0; i < 3; ++i)
            intersection[i] = v0[i] + u * e1[i] + v * e2[i];

        if ((u == 0. && v == 0.) || (u == 0. && v == 1.) || (u == 1. && v == 0.))
            // vertex
            return 1. / 6.;
        else if (u == 0. || v == 0. || u + v == 1.)
            // edge
            return 0.5;
        else
            // interior
            return 1.;
    }
}

} // namespace Bempp
