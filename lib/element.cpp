#include "element.hpp"

void TriangularElement::jacobian(
  int n, const Point* pts, double* result) const
{
  const Point* v = vertices;
  const double b[2][2] =
    {{2. * (-v[0].x + v[1].x),
      2. * (-v[1].x + v[2].x)},
     {2. * (-v[0].y + v[1].y),
      2. * (-v[1].y + v[2].y)}};

  const double J = 1/16. * (b[0][0] * b[1][1] - b[0][1] * b[1][0]);
  for (int i = 0; i < n; ++i)
    result[i] = J;
}

void QuadrilateralElement::jacobian(
  int n, const Point* pts, double* result) const
{
  // dx_i/du_j = 0.25 * (a_i * u_~j + b_{ij})

  const Point* v = vertices;
  const double a[2] = 
    {v[0].x - v[1].x + v[2].x - v[3].x,
     v[0].y - v[1].y + v[2].y - v[3].y};
  const double b[2][2] = 
    {{-v[0].x + v[1].x + v[2].x - v[3].x,
      -v[0].x - v[1].x + v[2].x + v[3].x},
     {-v[0].y + v[1].y + v[2].y - v[3].y,
      -v[0].y - v[1].y + v[2].y + v[3].y}};

  for (int i = 0; i < n; ++i)
    result[i] = 1/16. * ((a[0] * pts[i].y + b[0][0]) *
                         (a[1] * pts[i].x + b[1][1]) -
                         (a[0] * pts[i].x + b[0][1]) *
                         (a[1] * pts[i].y + b[1][0]));
}
