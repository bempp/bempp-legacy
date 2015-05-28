// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMMON_IMPL_HPP
#define HMAT_COMMON_IMPL_HPP

#include "common.hpp"
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polytope_distance_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/MP_Float.h>

namespace hmat {

inline IndexSetType fillIndexRange(std::size_t start, std::size_t stop) {
  IndexSetType result;
  result.reserve(stop - start);
  for (int i = start; i < stop; ++i)
    result.push_back(i);
  return result;
}

inline double obbDistance(const std::vector<Point>& points1,
    const std::vector<Point>& points2)
{
  typedef CGAL::Polytope_distance_d_traits_3<CgalKernel, CGAL::Gmpzf, double> Traits;
  typedef CGAL::Polytope_distance_d<Traits> Polytope_distance;

  Polytope_distance pd(begin(points1),end(points1),begin(points2),end(points2));

  return std::sqrt(CGAL::to_double (pd.squared_distance_numerator()) / CGAL::to_double (pd.squared_distance_denominator()));
  //return std::sqrt(pd.squared_distance());
}

inline double clusterDiameter(const std::vector<Point>& points)
{
  typedef CGAL::Min_sphere_of_spheres_d_traits_3<CgalKernel,double> Traits;
  typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
  typedef Traits::Sphere Sphere;

  Min_sphere ms;

  for (const auto& point : points)
    ms.insert(Sphere(point,0));

  return 2.*ms.radius();
}

}




#endif
